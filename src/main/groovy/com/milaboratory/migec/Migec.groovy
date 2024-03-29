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

import java.util.jar.JarFile

def version = (getClass().classLoader.findResource(JarFile.MANIFEST_NAME).text =~
        /Implementation-Version: (.+)/)[0][1]

def printHelp = {
    println " MiGEC Pipeline V$version "
    println ""
    println "[Main pipeline]"
    println "Run as \$java -jar migec-${version}.jar SCRIPT_NAME arguments"
    println ""
    println "where SCRIPT_NAME is one of the following:"
    println "Checkout"
    println "CheckoutBatch"
    println "Histogram"
    println "Assemble"
    println "AssembleBatch"
    println "CdrBlast"
    println "CdrBlastBatch"
    println "FilterCdrBlastResults"
    println "FilterCdrBlastResultsBatch"
    println "CreateCdrHypermGraph"
    println "Report"
    println ""
    println ""
    println "[Miscellaneous]"
    println "Executed via classpath"
    println "Run as \$java -cp migec-${version}.jar SCRIPT_NAME arguments"
    println ""
    println "where SCRIPT_NAME is one of the following:"
    println "com.milaboratory.migec.GroupByCdr"
    println "com.milaboratory.migec.SplitSpikeClonotypes"
    println "com.milaboratory.migec.BacktrackSequence"
    println "com.milaboratory.migec.HistQ"
    println "com.milaboratory.migec.UmiFrequencyTable"
}

def getScript = { String scriptName ->
    switch (scriptName.toUpperCase()) {
        case "CHECKOUT":
            return new Checkout()
        case "CHECKOUTBATCH":
            return new CheckoutBatch()
        case "HISTOGRAM":
            return new Histogram()
        case "ASSEMBLE":
            return new Assemble()
        case "ASSEMBLEBATCH":
            return new AssembleBatch()
        case "CDRBLAST":
            return new CdrBlast()
        case "CDRBLASTBATCH":
            return new CdrBlastBatch()
        case "FILTERCDRBLASTRESULTS":
            return new FilterCdrBlastResults()
        case "FILTERCDRBLASTRESULTSBATCH":
            return new FilterCdrBlastResultsBatch()
        case "CREATECDRHYPERMGRAPH":
            return new CreateCdrHypermGraph()
        case "REPORT":
            return new Report()
        case "-H":
        case "H":
        case "-HELP":
        case "HELP":
            printHelp()
            println ""
            System.exit(0)
            break

        default:
            printHelp()
            println ""
            println "Unknown MAIN PIPELINE script $scriptName"
            System.exit(2)
    }
}

if (args.length == 0)
    printHelp()
else {
    def script = getScript(args[0])
    try {
        Util.CMD_LINE = "migec-${version}.jar ${args.join(" ")}"
        Util.run(script, args.length > 1 ? args[1..-1].join(" ") : "")
    } catch (Exception e) {
        println "[ERROR] $e.message, see _migec_error.log for details"
        new File("_migec_error.log").withWriterAppend { writer ->
            writer.println("[${new Date()}]")
            writer.println("[Script]")
            writer.println(args[0])
            writer.println("[CommandLine]")
            writer.println("executing $Util.CMD_LINE")
            writer.println("[Message]")
            writer.println(e.message)
            writer.println("[StackTrace-Short]")
            writer.println(e.stackTrace.findAll { it.toString().contains("com.milaboratory.migec") }.join("\n"))
            writer.println("[StackTrace-Full]")
            e.printStackTrace(new PrintWriter(writer))
        }
        System.exit(1)
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
