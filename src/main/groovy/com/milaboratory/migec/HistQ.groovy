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

import java.util.zip.GZIPInputStream

if (args.length < 5) {
    println "[ERROR] Too few arguments provided"
    println "Usage: HistQ inputFastqFile output from(e.g. 0) to(e.g. read length, 100)"
    System.exit(2)
}

def infiles = args[0].split(","), ofile = args[1]
int from = Integer.parseInt(args[2]), to = Integer.parseInt(args[3])

double[][] qhist = new double[42][to - from]
int n = 0

infiles.each { infile ->
    def reader = new BufferedReader(new InputStreamReader(infile.endsWith(".gz") ?
            new GZIPInputStream(new FileInputStream(infile)) : new FileInputStream(infile)))

    while (reader.readLine() != null) {
        reader.readLine()
        reader.readLine()
        def q = reader.readLine().toCharArray()
        for (int i = from; i < to; i++) {
            qhist[((int) q[i] - 33)][i - from]++
        }
        if (n++ % 1000000 == 0)
            println "[${new Date()}] Processed $n reads"
    }
}

if (new File(ofile).parentFile)
    new File(ofile).parentFile.mkdirs()

new File(ofile).withPrintWriter { pw ->
    pw.println(qhist.collect { it.collect().join("\t") }.join("\n"))
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
