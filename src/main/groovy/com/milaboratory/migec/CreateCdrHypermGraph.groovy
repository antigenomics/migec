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

import java.util.concurrent.atomic.AtomicLong

import static com.milaboratory.migec.Util.BLANK_FIELD

def T = "10000"
def cli = new CliBuilder(usage: 'CreateCdrHypermGraph file_with_cdrs path/to/output/dir')
cli.h("usage")
cli.p(args: 1, 'Number of threads. Default: all available processors.')
cli.t(args: 1, "Clonotype size ratio threshold to record a hypermutation. Default: $T.")
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
def scriptName = getClass().canonicalName
def threshold = Math.abs(Math.log10(Double.parseDouble(opt.t ?: T)))
def inputFileName = opt.arguments()[0], outputPath = opt.arguments()[1]

int PERC_COL = 1, NT_SEQ_COL = 2, AA_SEQ_COL = 3, V_GENE_COL = 4, J_GENE_COL = 5,
    V_END_COL = 7, J_START_COL = 10

def cdr3Map = new HashMap()

println "[${new Date()} $scriptName] Reading in clonotypes"
new File(inputFileName).splitEachLine("\t") { line ->
    if (line[0].isInteger())
        cdr3Map.put(line[NT_SEQ_COL], [Double.parseDouble(line[PERC_COL]), line[AA_SEQ_COL],
                                       line[V_GENE_COL], line[J_GENE_COL],
                                       line[V_END_COL], line[J_START_COL]])
}

def hypermGraph = Collections.synchronizedMap(new HashMap())
def count = new AtomicLong()
println "[${new Date()} $scriptName] Searching one-mismatch pairs"
GParsPool.withPool THREADS, {
    cdr3Map.eachParallel { Map.Entry<String, Double> entry ->
        def thisSeq = entry.key
        def chars = thisSeq.toCharArray()
        def oldChar
        for (int i = 0; i < chars.length; i++) {
            oldChar = chars[i]
            // hash-based 1-loop single-mm search
            [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
                if (nt != oldChar) {
                    chars[i] = nt
                    def otherSeq = new String(chars)
                    def otherSeqEntry = cdr3Map[otherSeq]
                    if (otherSeqEntry != null &&  // mismatch exists
                            hypermGraph[otherSeq] != thisSeq) // not added already
                        if (Math.abs(Math.log10(otherSeqEntry[0] / entry.value[0])) < threshold)
                            hypermGraph.put(thisSeq, otherSeq)
                }
            }
            chars[i] = oldChar
        }
        def cc = count.incrementAndGet()
        if (cc % 50000 == 0)
            println "[${new Date()} $scriptName] $cc processed of ${cdr3Map.size()}"
    }
}
println "[${new Date()} $scriptName] Found ${hypermGraph.size()} hypermutations"

new File(outputPath).mkdirs()

println "[${new Date()} $scriptName] Writing node data"
new File(outputPath + "/nodes.txt").withPrintWriter { pw ->
    pw.println("NT_SEQ\tSHARE\tAA_SEQ\tV_SEGMENT\tJ_SEGMENT\tV_END\tJ_START")
    cdr3Map.each {
        pw.println(it.key + "\t" + Math.log10(it.value[0]) + "\t" + it.value[1..-1].collect().join("\t"))
    }
}

println "[${new Date()} $scriptName] Writing net and edge data"
new File(outputPath + "/net.txt").withPrintWriter { pw1 ->
    new File(outputPath + "/edges.txt").withPrintWriter { pw ->
        pw.println("SEQ1-SEQ2\tMISSENSE\tSUBSTITUTION\tSUBSTITUTION_CODON")
        hypermGraph.each {
            String fromNt = it.key, toNt = it.value
            pw1.println(fromNt + "\t" + toNt)
            String aaseq1 = cdr3Map[fromNt][1], aaseq2 = cdr3Map[toNt][1]
            def aaSubst = BLANK_FIELD, codonSubst = BLANK_FIELD
            for (int i = 0; i < aaseq1.size(); i++) {
                String from = aaseq1.charAt(i), to = aaseq2.charAt(i)
                if (from != to) {
                    aaSubst = from + "<>" + to
                }
            }
            for (int i = 0; i < fromNt.size(); i++) {
                String from = fromNt.charAt(i), to = toNt.charAt(i)
                if (from != to) {
                    int ntPos = i - i % 3
                    codonSubst = fromNt.substring(ntPos, ntPos + 3) + "<>" + toNt.substring(ntPos, ntPos + 3)
                }
            }
            pw.println(it.key + " (pp) " + it.value + "\t" + (aaSubst != BLANK_FIELD) + "\t" + aaSubst + "\t" + codonSubst)
        }
    }
}
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
