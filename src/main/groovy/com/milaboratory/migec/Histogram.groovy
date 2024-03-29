/*
 * Copyright (c) 2014-2024, OOO «MiLaboratory»
 *
 * IN NO EVENT SHALL THE INVENTORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE LICENSOR HAS NO
 * OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS. THE LICENSOR MAKES NO REPRESENTATIONS AND EXTENDS NO
 * WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
 * PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * TRADEMARK OR OTHER RIGHTS.
 */

package com.milaboratory.migec

import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray

def DEFAULT_PARENT_CHILD_RATIO = "1"

def cli = new CliBuilder(usage: 'Histogram [options] checkout_dir/ output_dir/')
cli.h("usage")
cli.q(args: 1, argName: 'read quality (phred)', "barcode region quality threshold. " +
        "Default: $Util.DEFAULT_UMI_QUAL_THRESHOLD")
cli.p(args: 1, 'number of threads to use')
cli._(longOpt: 'only-first-read',
        'Use only first read (as they were in raw FASTQ) ' +
                'which typically has higher quality to compute statistics.')
cli._(longOpt: "collision-ratio", args: 1, argName: "double, <= 1.0",
        "Min parent-to-child MIG size ratio to count collision events. Default value: $DEFAULT_PARENT_CHILD_RATIO")

def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
byte umiQualThreshold = opt.q ? Byte.parseByte(opt.q) : Util.DEFAULT_UMI_QUAL_THRESHOLD
def inputDir = opt.arguments()[0]
boolean onlyFirstRead = opt.'only-first-read'
double collisionRatioThreshold = Double.parseDouble(opt.'collision-ratio' ?: DEFAULT_PARENT_CHILD_RATIO)

def scriptName = getClass().canonicalName
int nBins = 17

def scale = { Integer value ->
    (int) (Math.min(Math.log((double) value) / Math.log(2.0), nBins - 1))
}

def samples = new File("$inputDir/checkout.filelist.txt").readLines().findAll { !it.startsWith("#") }.collect {
    it.split("\t")
}

def outputDir = opt.arguments()[1]
new File(outputDir).mkdirs()

def bins = (0..<nBins).collect { (int) Math.pow(2, it) }

def BASE_HEADER = "#SAMPLE_ID\tSAMPLE_TYPE",
    HEADER = BASE_HEADER + "\t" + bins.join("\t"),
    ESTIMATES_HEADER = BASE_HEADER + "\t" +
            "TOTAL_READS\tTOTAL_MIGS\t" +
            "OVERSEQ_THRESHOLD\tCOLLISION_THRESHOLD\t" +
            "UMI_QUAL_THRESHOLD\tUMI_LEN"

def maxUmiSz = -1
def umiPwmUnitsTotal = new AtomicIntegerArray(4000), umiPwmReadsTotal = new AtomicIntegerArray(4000),
    umiGoodUnitsTotal = new AtomicInteger(), umiGoodReadsTotal = new AtomicInteger()

new File("$outputDir/overseq.txt").withPrintWriter { oWriter ->
    oWriter.println(HEADER)

    new File("$outputDir/overseq-units.txt").withPrintWriter { ouWriter ->
        ouWriter.println(HEADER)

        new File("$outputDir/collision1.txt").withPrintWriter { cWriter ->
            cWriter.println(HEADER)

            new File("$outputDir/collision1-units.txt").withPrintWriter { cuWriter ->
                cuWriter.println(HEADER)


                new File("$outputDir/estimates.txt").withPrintWriter { eWriter ->
                    eWriter.println(ESTIMATES_HEADER)

                    new File("$outputDir/pwm-units.txt").withPrintWriter { pwPwmUnits ->

                        new File("$outputDir/pwm.txt").withPrintWriter { pwPwm ->

                            samples.each { sampleEntry ->
                                println "[${new Date()} $scriptName] Processing ${sampleEntry[0]} (${sampleEntry[1]})"

                                def umiPwmUnits = new AtomicIntegerArray(4000), umiPwmReads = new AtomicIntegerArray(4000),
                                    umiGoodUnits = new AtomicInteger(), umiGoodReads = new AtomicInteger()
                                int umiSz = -1

                                // Accumulate UMIs
                                def umiCountMap = new HashMap<String, Integer>()
                                def reader = Util.getReader(sampleEntry[2])

                                String header

                                int nReads = 0
                                while ((header = reader.readLine()) != null) {
                                    def umi = Util.getUmi(header, umiQualThreshold)

                                    if (umi != null && (!onlyFirstRead || header.contains(" R1 "))) {
                                        umiCountMap.put(umi, (umiCountMap.get(umi) ?: 0) + 1)

                                        if (umiSz < 0)
                                            umiSz = umi.length()
                                        else if (umiSz != umi.length()) {
                                            println "[ERROR] UMIs of various sizes are not supported in the same sample"
                                            System.exit(2)
                                        }
                                    }

                                    reader.readLine()
                                    reader.readLine()
                                    reader.readLine()

                                    if (++nReads % 1000000 == 0)
                                        println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs so far"
                                }

                                maxUmiSz = Math.max(maxUmiSz, umiSz)

                                println "[${new Date()} $scriptName] Processed $nReads, ${umiCountMap.size()} UMIs total"

                                def overseqHist = new AtomicIntegerArray(nBins), overseqHistUnits = new AtomicIntegerArray(nBins),
                                    collisionHist = new AtomicIntegerArray(nBins), collisionHistUnits = new AtomicIntegerArray(nBins)

                                def nUmis = new AtomicInteger()

                                GParsPool.withPool THREADS, {
                                    umiCountMap.eachParallel { Map.Entry<String, Integer> umiEntry ->
                                        int thisCount = umiEntry.value, bin = scale(thisCount)

                                        // Append to cumulative overseq
                                        overseqHist.addAndGet(bin, thisCount)
                                        overseqHistUnits.incrementAndGet(bin)

                                        // Calculate 1-mm collisions
                                        char[] umi = umiEntry.key.toCharArray()
                                        boolean good = true

                                        //loops:
                                        for (int i = 0; i < umi.length; i++) {
                                            char prevChar = umi[i]
                                            for (int j = 0; j < 4; j++) {
                                                char nt = Util.NTS[j]
                                                if (prevChar != nt) {
                                                    umi[i] = nt
                                                    def otherCount = umiCountMap.get(new String(umi))
                                                    if (otherCount != null &&
                                                            thisCount < collisionRatioThreshold * otherCount) {
                                                        collisionHist.addAndGet(bin, thisCount)
                                                        collisionHistUnits.incrementAndGet(bin)
                                                        good = false
                                                        break
                                                        //break loops, let's not over-count collisions
                                                    }
                                                }
                                            }

                                            // Looks like groovy compiler sometimes fails
                                            // with break-label
                                            if (!good)
                                                break
                                            
                                            umi[i] = prevChar
                                        }

                                        if (good) {
                                            umiGoodReads.addAndGet(thisCount)
                                            umiGoodUnits.incrementAndGet()
                                            umiGoodReadsTotal.addAndGet(thisCount)
                                            umiGoodUnitsTotal.incrementAndGet()
                                            for (int i = 0; i < umi.length; i++) {
                                                int x = 4 * i + Util.nt2code(umi[i])
                                                umiPwmUnits.incrementAndGet(x)
                                                umiPwmReads.addAndGet(x, thisCount)
                                                umiPwmUnitsTotal.incrementAndGet(x)
                                                umiPwmReadsTotal.addAndGet(x, thisCount)
                                            }
                                        }

                                        int nUmisCurrent = nUmis.incrementAndGet()

                                        if (nUmisCurrent % 10000 == 0)
                                            println "[${new Date()} $scriptName] Collecting stats, $nUmisCurrent UMIs processed"
                                    }
                                }

                                def row = sampleEntry[0..1].join("\t")

                                def estimateHist = new double[nBins]

                                // check if there is any over-sequencing at all
                                boolean overSeqPresent = (0..2).sum { overseqHist.get(it) } <
                                        2.0 * (3..<nBins).sum { overseqHist.get(it) }

                                int collThreshold, overseqThreshold

                                if (overSeqPresent) {
                                    // Find peak
                                    // smooth first
                                    (0..<nBins).each { i ->
                                        def avg = 0, n = 0
                                        (-1..1).each { j ->
                                            int ii = i + j
                                            if (ii > -1 && ii < nBins) {
                                                avg += overseqHist.get(ii)
                                                n++
                                            }
                                        }
                                        estimateHist[i] = avg / n
                                    }
                                    // look for maximum in appropriate range
                                    def overseqPeak = (3..10).max { estimateHist[it] }
                                    overseqThreshold = Math.round(Math.pow(2.0, overseqPeak / 2.0))
                                } else {
                                    overseqThreshold = 1 // assemble all
                                }

                                collThreshold = overseqThreshold // for now

                                oWriter.println(row + "\t" + Util.toString(overseqHist))
                                cWriter.println(row + "\t" + Util.toString(collisionHist))
                                ouWriter.println(row + "\t" + Util.toString(overseqHistUnits))
                                cuWriter.println(row + "\t" + Util.toString(collisionHistUnits))
                                eWriter.println(row + "\t" +
                                        nReads + "\t" + nUmis + "\t" +
                                        overseqThreshold + "\t" + collThreshold + "\t" +
                                        umiQualThreshold + "\t" + umiSz)


                                def PWM_HEADER = ESTIMATES_HEADER = BASE_HEADER + "\tNT\t" + (1..umiSz).join("\t")
                                pwPwmUnits.println(PWM_HEADER)
                                pwPwm.println(PWM_HEADER)
                                for (int i = 0; i < 4; i++) {
                                    pwPwmUnits.print(row + "\t" + Util.code2nt(i))
                                    pwPwm.print(row + "\t" + Util.code2nt(i))
                                    for (int j = 0; j < umiSz; j++) {
                                        int x = 4 * j + i
                                        pwPwmUnits.print("\t" + umiPwmUnits.get(x) / (double) umiGoodUnits.get())
                                        pwPwm.print("\t" + umiPwmReads.get(x) / (double) umiGoodReads.get())
                                    }
                                    pwPwmUnits.println()
                                    pwPwm.println()
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

new File("$outputDir/pwm-summary-units.txt").withPrintWriter { pwPwmUnits ->

    new File("$outputDir/pwm-summary.txt").withPrintWriter { pwPwm ->
        def PWM_HEADER = ESTIMATES_HEADER = "#NT\t" + (1..maxUmiSz).join("\t")
        pwPwmUnits.println(PWM_HEADER)
        pwPwm.println(PWM_HEADER)
        for (int i = 0; i < 4; i++) {
            pwPwmUnits.print(Util.code2nt(i))
            pwPwm.print(Util.code2nt(i))
            for (int j = 0; j < maxUmiSz; j++) {
                int x = 4 * j + i
                pwPwmUnits.print("\t" + umiPwmUnitsTotal.get(x) / (double) umiGoodUnitsTotal.get())
                pwPwm.print("\t" + umiPwmReadsTotal.get(x) / (double) umiGoodReadsTotal.get())
            }
            pwPwmUnits.println()
            pwPwm.println()
        }
    }
}

Util.printCmd("$outputDir/histogram.cmd.txt")
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
NA
