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

import static com.milaboratory.migec.Util.BLANK_PATH

def cli = new CliBuilder(usage:
        "CheckoutBatch [options, see Checkout] barcode_file output_dir/")
cli.h("usage")

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

def options = args.size() > 2 ? args[0..-3] : [], barcodesFileName = args[-2], outputDir = args[-1]

def scriptName = getClass().canonicalName

if (!options.contains("--append"))
    options.add("--append")

boolean anyMissing = false
def filesHash = new HashSet<String>()
new File(barcodesFileName).splitEachLine("\t") { splitLine ->
    if (splitLine.size() > 0 &&
            !splitLine[0].startsWith("#") && splitLine.size() > 3 && splitLine[3] != BLANK_PATH) {
        def fastq1 = splitLine[3], fastq2 = splitLine.size() > 4 ? splitLine[4] : BLANK_PATH
        filesHash.add([fastq1, fastq2].join(" "))

        [new File(fastq1), new File(fastq2)].each {
            if (!it.exists()) {
                println "${it.absolutePath} not found"
                anyMissing = true
            }
        }
    }
}

if (anyMissing) {
    println "[ERROR] There were some missing FASTQ files.. Terminating"
    System.exit(2)
}

def fileList = filesHash.collect().sort()

println "[${new Date()} $scriptName] Will run Checkout for the following fastq pairs:\n" +
        "${fileList.join("\n")}\n" +
        "barcodes: $barcodesFileName\n" +
        "output path: $outputDir"

new File(outputDir).deleteDir()

// Clear existing logs
[new File("$outputDir/checkout.filelist.txt"), new File("$outputDir/checkout.log.txt")].each {
    if (it.exists())
        it.delete()
}

new File(outputDir).mkdirs()

fileList.each {
    Util.run(new Checkout(), [options, barcodesFileName, it, outputDir].flatten().join(" "))
}

Util.printCmd(outputDir + "/checkout.cmd.txt")
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
