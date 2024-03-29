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

import java.nio.file.FileSystems
import java.nio.file.Files

def homeDir = "."

new File("$homeDir/summary/").deleteDir()
new File("$homeDir/summary/").mkdir()

boolean firstSample = true
new File(homeDir).eachDir { sampleDir ->
    def sampleName = sampleDir.name

    println "Processing $sampleName"

    sampleDir.eachDir { outputDir ->
        def outputName = outputDir.name

        // copy logs

        def logFile = outputDir.listFiles().find { it.name.endsWith(".log.txt") }

        if (logFile != null) {
            new File("$homeDir/summary/${outputName}.log.txt").withWriterAppend { w ->
                def logLines = logFile.readLines()

                if (firstSample)
                    w.println(logLines[0])

                logLines[1..-1].each { w.println(it) }
            }
        }

        // copy final clonotype output

        if (outputName == "cdrfinal")
            Files.copy(
                    FileSystems.default.getPath(outputDir.listFiles().find {
                        it.name.endsWith(".filtered.cdrblast.txt")
                    }.absolutePath),
                    FileSystems.default.getPath("$homeDir/summary/${sampleName}.filtered.cdrblast.txt")
            )

        // histograms - just stack all of them
        if (outputName == "histogram") {
            outputDir.listFiles().each { histFile ->
                new File("$homeDir/summary/$histFile.name").withWriterAppend { w ->
                    def histLines = histFile.readLines()

                    if (firstSample)
                        w.println(histLines[0])

                    histLines[1..-1].each { w.println(it) }
                }
            }
        }
    }

    firstSample = false
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
