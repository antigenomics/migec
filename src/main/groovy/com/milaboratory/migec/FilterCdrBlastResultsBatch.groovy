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

def R_A_T = "1.0", S_F_R = "20.0"
def scriptName = getClass().canonicalName
def cli = new CliBuilder(usage: "$scriptName [options] cdrblast_dir/ output_dir/")
cli.h("usage")
cli.r(args: 1, argName: 'read accumulation threshold',
        "Only clonotypes that have a ratio of (reads after correction) / " +
        "(uncorrected reads) greater than that threshold are retained. Default: $R_A_T")
cli._(longOpt: 'collapse', "Collapse clonotypes by CDR3 and use top V and J chains")
cli.p(args: 1, "number of threads to use. Default: all available processors")
cli.s(longOpt: "singleton-filter",
        "Perform frequency-based filtering of singletons, i.e. clonotypes that are represented by a single MIG")
cli._(longOpt: "singleton-filter-ratio", args: 1, argName: "float>1.0",
        "Parent-to-child ratio for frequency-based filtering of singleton clonotypes [default = $S_F_R]")
cli.n("Include non-functional CDR3s")
cli.c("Include CDR3s that do not begin with a conserved C or end with a conserved W/F")

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

def baseArgs = [["-r", opt.r ?: R_A_T]]
if (opt.s)
    baseArgs = [baseArgs, ["-s"]]
if (opt.'--singleton-filter-ratio')
    baseArgs = [baseArgs, ["--singleton-filter-ratio", opt.'--singleton-filter-ratio' ?: S_F_R]]
if (opt.n)
    baseArgs = [baseArgs, ["-n"]]
if (opt.c)
    baseArgs = [baseArgs, ["-c"]]
if (opt.p)
    baseArgs = [baseArgs, ["-p", opt.p.toString()]]
if (opt.'collapse')
    baseArgs = [baseArgs, ["--collapse"]]

String inputDir = opt.arguments()[0], outputDir = opt.arguments()[1]

def cdrBlastBatchFile = new File("$inputDir/cdrblast.log.txt")
if (!cdrBlastBatchFile.exists()) {
    println "[ERROR] CdrBlast log file (cdrblast.log.txt) not found in input dir"
    System.exit(2)
}

new File(outputDir).mkdirs()

// collect raw-assembled pairs
def rawSampleMap = new HashMap(), assembledSampleMap = new HashMap()

cdrBlastBatchFile.splitEachLine("\t") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        if (splitLine[1].toLowerCase() == "asm") {
            assembledSampleMap.put(splitLine[0], splitLine[2])
        } else if (splitLine[1].toLowerCase() == "raw") {
            rawSampleMap.put(splitLine[0], splitLine[2])
        } else {
            println "[ERROR] Unknown sample type ${splitLine[1]}. Either \'raw\' or \'asm\' expected"
            System.exit(2)
        }

        if (!new File(splitLine[2]).exists()) {
            println "[ERROR] File not found for sample ${splitLine[0..2]}."
            System.exit(2)
        }
    }
}

def onlyInRaw = rawSampleMap.keySet().findAll { !assembledSampleMap[it] },
    onlyInAssembled = assembledSampleMap.keySet().findAll { !rawSampleMap[it] }

if (onlyInRaw.size() > 0) {
    println "[WARNING] The following samples have only raw data analyzed by CdrBlast and will be skipped: " +
            "${onlyInRaw.join(", ")}"
}
if (onlyInAssembled.size() > 0) {
    println "[WARNING] The following samples have only assembled data analyzed by CdrBlast and will be skipped: " +
            "${onlyInAssembled.join(", ")}"
}

if (onlyInRaw.size() == rawSampleMap.size()) {
    println "[ERROR] No samples have both raw and assembled data analyzed by CdrBlast. Terminating"
    System.exit(2)
}

println "[${new Date()} $scriptName] Running batch hot-spot error filtering for CdrBlast results.."

// START LOGGING
def logFile = new File(outputDir + "/cdrblastfilter.log.txt")
if (logFile.exists())
    logFile.delete()
else
    logFile.absoluteFile.parentFile.mkdirs()

logFile.withPrintWriter { pw ->

    pw.println(Util.CDRBLASTFILTER_LOG_HEADER)

    assembledSampleMap.each {
        def sampleId = it.key, asmFileName = it.value, rawFileName = rawSampleMap[sampleId]

        if (rawFileName) {
            def finalArgs = [baseArgs, [asmFileName, rawFileName, "$outputDir/${sampleId}.filtered.cdrblast.txt"]]
            def logLine = Util.run(new FilterCdrBlastResults(), finalArgs.flatten().join(" "))
            pw.println("$sampleId\t" + logLine)
        }
    }
}

Util.printCmd(outputDir + "/cdrblastfilter.cmd.txt")
if (opt.n)
    new File(outputDir + "/cdrblastfilter.nc.txt").withPrintWriter { it.println() }
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
