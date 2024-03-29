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

def cli = new CliBuilder(usage: 'SplitSpikeClonotypes file_with_spikes file_with_cdrs output_prefix')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}
def spikeFileName = opt.arguments()[0], inputFileName = opt.arguments()[1], outputFilePrefix = opt.arguments()[2]

def spikeList = new File(spikeFileName).readLines()
def spikeHash = new HashSet<String>(spikeList)
// add 1 mm variants also
spikeList.each {
    def chars = it.toCharArray()
    def oldChar
    for (int i = 0; i < chars.length; i++) {
        oldChar = chars[i]
        [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
            if (nt != oldChar) {
                chars[i] = nt
                spikeHash.add(new String(chars))
            }
        }
        chars[i] = oldChar
    }
}

// Iterate through file
int NT_SEQ_COL = 2
new File(outputFilePrefix + ".sample.txt").withPrintWriter { pw1 ->
    new File(outputFilePrefix + ".spikes.txt").withPrintWriter { pw2 ->
        new File(inputFileName).splitEachLine("\t") { line ->
            if (line[0].isInteger()) {
                String seq = line[NT_SEQ_COL]
                def chars = seq.toCharArray()
                def oldChar
                boolean isSpike = false
                for (int i = 0; i < chars.length; i++) {
                    oldChar = chars[i]
                    // two mismatches total
                    [(char) 'A', (char) 'T', (char) 'G', (char) 'C'].each { char nt ->
                        chars[i] = nt
                        if (spikeHash.contains(new String(chars)))
                            isSpike = true
                    }
                    chars[i] = oldChar
                    if (isSpike)
                        break
                }
                if (isSpike)
                    pw2.println(line.join("\t"))
                else
                    pw1.println(line.join("\t"))
            } else {
                pw1.println(line.join("\t"))
                pw2.println(line.join("\t"))
            }
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
