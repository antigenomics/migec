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

//
// Internal routine checks how many different CDR3 sequences are there per a single MIG
// To be applied to the output of CdrBlast processing of raw data with --cdr3-fastq-file option
// CDR3s with a single mismatch are collapsed, as here we're interested only in spotting events when
// two different clonotypes got marked by the same UMI, not errors
//

package com.milaboratory.migec

def (inputFileName, outputFileName) = args[0..1]

def cdrCountsByUmi = new HashMap<String, HashMap<String, Integer>>()

println "[${new Date()}] Started reading"
new File(inputFileName).withReader { reader ->
    def header
    def counter = 0

    while ((header = reader.readLine()) != null) {
        def splitHeader = header.split(" "),
            umi = splitHeader.find { it.startsWith("UMI:") }.split(":")[1],
            cdr = splitHeader.find { it.startsWith("CDR3:") }.split(":")[1]

        def cdrCounts = cdrCountsByUmi[umi]
        if (!cdrCounts)
            cdrCountsByUmi.put(umi, cdrCounts = new HashMap<String, Integer>())
        cdrCounts.put(cdr, (cdrCounts[cdr] ?: 0) + 1)

        reader.readLine()
        reader.readLine()
        reader.readLine()

        if (++counter % 100000 == 0)
            println "[${new Date()}] $counter reads processed"
    }
}

println "[${new Date()}] Processing and writing output"
new File(outputFileName).withPrintWriter { pw ->
    pw.println("#umi\ttotal\tcdr3\tcount")

    int counter = 0

    cdrCountsByUmi.sort { -it.value.values().sum() }.each {
        def umi = it.key, cdrCounts = it.value
        int total = cdrCounts.values().sum()

        def cdrCountsCollapsed = new HashMap<String, Integer>()

        // collapse by 1mm
        cdrCounts.each {
            def cdr = it.key, cdrChars = cdr.toCharArray(),
                count = it.value
            boolean bad = false

            for (int i = 0; i < cdrChars.length; i++) {
                char oldChar = cdrChars[i]
                [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
                    if (nt != oldChar && !bad) {
                        cdrChars[i] = nt
                        def otherCdr = new String(cdrChars)
                        def otherCount = cdrCounts[otherCdr]

                        if (count < otherCount) {
                            cdrCountsCollapsed.put(otherCdr, (cdrCountsCollapsed[otherCdr] ?: 0) + count)
                            bad = true
                        }
                    }
                }
                cdrChars[i] = oldChar
                if (bad)
                    break
            }

            if (!bad)
                cdrCountsCollapsed.put(cdr, (cdrCountsCollapsed[cdr] ?: 0) + count)
        }

        cdrCountsCollapsed.sort { it.key }.sort { -it.value }.each {
            pw.println([umi, total, it.key, it.value].join("\t"))
        }

        if (++counter % 1000 == 0)
            println "[${new Date()}] $counter MIGs processed"
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
