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

import com.milaboratory.migec.alignment.ReadData

def DEFAULT_THRESHOLD = "10"

def cli = new CliBuilder(usage:
        'FilterReadsByCount [options] input.fastq[.gz] output.fastq[.gz]\n[for benchmarking purposes]')
cli.h("usage")
cli.t(args: 1, "Count threshold. $DEFAULT_THRESHOLD")
cli.g("Grouped output: reads collapsed, quality average, count in header.")

def scriptName = getClass().canonicalName

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

int threshold = (opt.t ?: DEFAULT_THRESHOLD).toInteger()
boolean group = opt.g ? true : false
def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

def reader = Util.getReader(inputFileName)


def readDataMap = new HashMap<String, ReadData>()

int nReads = 0

while (reader.readLine() != null) {
    def seq = reader.readLine()
    reader.readLine()
    def qual = reader.readLine()

    def readData = readDataMap[seq]

    if (readData == null)
        readDataMap.put(seq, readData = new ReadData(seq))

    readData.append(qual, group)

    if (++nReads % 500000 == 0)
        println "[${new Date()} $scriptName] Loaded $nReads reads"
}

nReads = 0
int nPassedReads = 0

def writer = Util.getWriter(outputFileName)

if (group) {
    int nReadGroups = 0
    readDataMap.each {
        int count = it.value.count
        if (count >= threshold) {
            writer.writeLine("@GroupedRead Count:$it.value.count\n" + it.key + "\n+\n" + it.value.finalizeQual())
            nPassedReads += count
        }

        nReads += it.value.count

        if (++nReadGroups % (500000 / threshold) == 0)
            println "[${new Date()} $scriptName] Scanned $nReads reads, $nPassedReads reads passed"
    }
} else {
    reader = Util.getReader(inputFileName)

    def header

    while ((header = reader.readLine()) != null) {
        def seq = reader.readLine()

        int count = readDataMap[seq].count

        if (count >= threshold) {
            reader.readLine()
            def qual = reader.readLine()
            writer.writeLine(header + "\n" + seq + "\n+\n" + qual)
            nPassedReads++
        } else {
            reader.readLine()
            reader.readLine()
        }

        if (++nReads % 500000 == 0)
            println "[${new Date()} $scriptName] Scanned $nReads reads, $nPassedReads reads passed"
    }
}

writer.close()
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
