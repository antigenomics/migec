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

import static com.milaboratory.migec.Util.BLANK_FIELD

def DEFAULT_ASSEMBLE_MASK = "1:1", DEFAULT_PARENT_CHILD_RATIO = "0.1",
    DEFAULT_ASSEMBLY_OFFSET = "5", DEFAULT_ASSEMBLY_MISMATCHES = "5", DEFAULT_ASSEMBLY_ANCHOR = "10"
def cli = new CliBuilder(usage: 'AssembleBatch [options] checkout_dir/ histogram_dir/ output_dir/')
cli.h("usage")
cli.p(args: 1, 'number of threads to use')
cli._(longOpt: 'default-mask', args: 1, "Mask, default for all samples, see --sample-metadata. " +
        "Default: $DEFAULT_ASSEMBLE_MASK.")
cli._(longOpt: 'sample-metadata', args: 1, argName: 'file name',
        "A tab-delimited file indicating which samples to process and containing three columns:\n" +
                "sample_name (tab) file_type (tab) mask\n" +
                "Allowed file types\n" +
                "paired, unpaired, overlapped\n" +
                "Mask column is optional and will be ignored for unpaired and overlapped reads. " +
                "Mask indicates read(s) in pair that should be assembled. " +
                "Allowed values are\n1:0 (R1 assembled), 0:1 (R2 assembled) and 1:1 (both reads assembled, [default])")
cli._(longOpt: 'force-collision-filter', "Forced collision filtering for all samples.")
cli._(longOpt: 'force-overseq', args: 1, argName: 'int', "Forced minimal MIG size for all samples.")
cli._(longOpt: "collision-ratio", args: 1, argName: "double, < 1.0",
        "Min parent-to-child MIG size ratio for collision filtering. Default value: $DEFAULT_PARENT_CHILD_RATIO")
cli._(longOpt: "assembly-offset", args: 1, argName: "integer",
        "Assembly offset range. Default: $DEFAULT_ASSEMBLY_OFFSET")
cli._(longOpt: "assembly-mismatches", args: 1, argName: "integer",
        "Assembly max mismatches. Default: $DEFAULT_ASSEMBLY_MISMATCHES")
cli._(longOpt: "assembly-anchor", args: 1, argName: "integer",
        "Assembly anchor region half size. Default: $DEFAULT_ASSEMBLY_ANCHOR")
cli._(longOpt: 'only-first-read',
        'Use only first read (as they were in raw FASTQ), ' +
                'can improve assembly quality for non-oriented reads when' +
                'second read quality is very poor.')
cli.c("compressed output")

def scriptName = getClass().canonicalName
def opt = cli.parse(args)

if (opt == null || opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}
if (opt.h) {
    cli.usage()
    System.exit(0)
}

def checkoutDir = opt.arguments()[0], histogramDir = opt.arguments()[1],
    outputPath = opt.arguments()[2]

def defaultMask = opt.'default-mask' ?: DEFAULT_ASSEMBLE_MASK,
    forcedCollFilter = opt.'force-collision-filter',
    forcedOverseq = opt.'force-overseq' ? (opt.'force-overseq').toString().toInteger() : null

def checkoutSampleListFile = new File(checkoutDir + "/checkout.filelist.txt")
if (!checkoutSampleListFile.exists()) {
    println "[ERROR] Sample list (checkout.filelist.txt) absent in Checkout folder $checkoutDir"
    System.exit(2)
}
def histogramEstimatesFile = new File(histogramDir + "/estimates.txt")
if (!histogramEstimatesFile.exists()) {
    println "[ERROR] Histogram output (estimates.txt) absent in Histogram folder $histogramDir"
    System.exit(2)
}

new File(outputPath).mkdirs()

// Read filter
Map sampleFilter = null
String sampleFilterFileName = opt.'sample-metadata' ?: null

// File names from Checkout output
def sampleFileNamesMap = checkoutSampleListFile.readLines().findAll {
    !it.startsWith("#")   // no headers
}.collectEntries {
    def splitLine = it.split('\t')
    def sampleName = splitLine[0], sampleType = splitLine[1]
    [("$sampleName\t$sampleType".toString()): splitLine[2..3]]
}

if (sampleFilterFileName) {
    sampleFilter = new HashMap<String, List<Integer>>()

    new File(sampleFilterFileName).splitEachLine("\t") { splitLine ->
        if (!splitLine[0].startsWith("#")) {
            def sampleName = splitLine[0], sampleType = splitLine[1].toLowerCase(),
                sampleKey = "$sampleName\t$sampleType".toString()

            if (!Util.FILE_TYPES.any { sampleType == it }) {
                println "[ERROR] Bad sample type $sampleType"
                System.exit(2)
            }

            if (!sampleFileNamesMap.containsKey(sampleKey)) {
                println "[ERROR] Sample list is inconsistent between Checkout log and filter: " +
                        "$sampleKey not present in Checkout log"
                System.exit(2)
            }

            def sampleMask = splitLine.size() > 2 ? splitLine[2] : defaultMask
            sampleMask = sampleMask == BLANK_FIELD ? defaultMask : sampleMask

            if (sampleType.toUpperCase() == "paired" && !Util.MASKS.contains(sampleMask)) {
                println "[ERROR] Bad assembly mask $sampleMask for paired reads. " +
                        "Allowed values: ${Util.MASKS.collect().join(", ")}"
            }

            sampleFilter.put("$sampleName\t$sampleType".toString(), sampleMask)
        }
    }
}

def baseArgs = []
if (opt.c)
    baseArgs.add(['-c'])
if (opt.'only-first-read')
    baseArgs.add(['--only-first-read'])
if (opt.p)
    baseArgs.add(['-p', opt.p])
if (opt.'collision-ratio')
    baseArgs.add(['--collision-ratio', opt.'collision-ratio'])
if (opt.'assembly-offset')
    baseArgs.add(['--assembly-offset', opt.'assembly-offset'])
if (opt.'assembly-mismatches')
    baseArgs.add(['--assembly-mismatches', opt.'assembly-mismatches'])
if (opt.'assembly-anchor')
    baseArgs.add(['--assembly-anchor', opt.'assembly-anchor'])

double collisionFactorThreshold = 0.05

def logFile = new File(outputPath + "/assemble.log.txt")
if (logFile.exists())
    logFile.delete()
else
    logFile.absoluteFile.parentFile.mkdirs()

logFile.withPrintWriter { pw ->
    pw.println(Util.ASSEMBLE_LOG_HEADER)

    println "[${new Date()} $scriptName] Starting batch assembly.."
    histogramEstimatesFile.readLines().findAll { !it.startsWith("#") }.each {   // skip header
        def splitLine = it.split('\t')
        def sampleName = splitLine[0], sampleType = splitLine[1], sampleKey = "$sampleName\t$sampleType".toString()

        if (!sampleFilter || sampleFilter.containsKey(sampleKey)) {
            def mask = sampleFilter ? sampleFilter[sampleKey] : defaultMask

            if (sampleType != "paired" || mask != "0:0") {
                // Parse threshold estimates, check if it is safe to filter collisions
                def totalUmis = splitLine[3].toInteger(), overseqThreshold = forcedOverseq ?: splitLine[4].toInteger(),
                    collThreshold = splitLine[5].toInteger(), umiQualThreshold = Byte.parseByte(splitLine[6]),
                    umiLen = splitLine[7].toInteger(), filterCollisions = false

                if (forcedCollFilter || (collThreshold >= overseqThreshold &&
                        totalUmis < collisionFactorThreshold * Math.pow(4, umiLen - 1))) // # collisions << # starting molecules
                    filterCollisions = true

                // More messy argument passing
                def assembleArgs = [
                        baseArgs,
                        ['-m', overseqThreshold],
                        ['-q', umiQualThreshold],
                        ['--mask', mask]
                ]

                if (filterCollisions)
                    assembleArgs.add(['--filter-collisions'])

                // Pass filenames for I/O
                def sampleFileNames = sampleFileNamesMap[sampleKey]

                if (!sampleFileNames) {
                    println "[ERROR] Sample list is inconsistent between Checkout log and Histogram output: " +
                            "$sampleKey not present in Checkout log"
                    System.exit(2)
                }

                assembleArgs.add(sampleFileNames)
                assembleArgs.add([outputPath])

                String stats = Util.run(new Assemble(), assembleArgs.flatten().join(" "))

                pw.println(sampleKey + "\t" + stats)
            }
        }
    }
}

Util.printCmd(outputPath + "/assemble.cmd.txt")
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
