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
import static com.milaboratory.migec.Util.BLANK_PATH

def scriptName = getClass().canonicalName

def cli = new CliBuilder(usage: "$scriptName [options] " +
        "[checkout_dir/ or ${BLANK_PATH}] [assemble_dir/ or ${BLANK_PATH}] output_dir/\n" +
        "Either --sample-metadata or -R argument is required.")
cli.h("usage")
cli.p(args: 1, "Number of threads to use")
cli._(longOpt: "sample-metadata", args: 1, argName: "file",
        "An optional tab-delimited file indicating which samples to process and containing 6 columns:\n" +
                "sample_id|species|gene|file_types|mask|quality_threshold_pair\n" +
                "Samples not included will be omitted by CdrBlast. " +
                "File types, mask and quality threshold columns could be either omitted or set to " - "")
cli._(args: 1, longOpt: "blast-path",
        "Path to blast executable.")
cli._(longOpt: "no-sort",
        "Do not sort output files, which will speed up processing. " +
                "Could be used for full pipeline as FilterCdrBlastResults will provide final clonotype table in sorted format.")
cli._(longOpt: "all-segments",
        "Use full V/D/J segment library (including pseudogens, etc).")
cli._(longOpt: "all-alleles",
        "Use full list of alleles (uses only major *01 alleles by default)")
cli._(longOpt: "print-library",
        "Prints out allowed species-gene pairs. " +
                "To account non-functional segment data use together with --all-segments")
cli._(longOpt: "default-mask", args: 1, argName: "R1=0/1:R2=0/1",
        "Mask, default for all samples, will be ignored for unpaired and overlapped reads, " +
                "allowed values are: 1:0 (R1 assembled), 0:1 (R2 assembled) and 1:1 (both reads assembled, default)")
cli.R(longOpt: "default-gene", args: 1, argName: "gene",
        "Gene, default for all samples, allowed values: 'TRA', 'TRB', 'TRG', 'TRD', 'IGL', 'IGK' or 'IGH'. " +
                "Use --print-library for the list of allowed species-gene combinations.")
cli.S(longOpt: "default-species", args: 1, argName: "species",
        "Species, default for all samples, allowed values: 'HomoSapiens'[default], 'MusMusculus', ... " +
                "Use --print-library for the list of allowed species-gene combinations.")
cli._(longOpt: "default-file-types", args: 1, argName: "type1,..",
        "Accepted file types, default for all samples, " +
                "comma-separated list of file types to be processed for a given sample, allowed values: ${Util.FILE_TYPES.join(", ")}")
cli.q(longOpt: "default-quality-threshold", args: 1, argName: "Phred,CQS",
        "Quality threshold pair, default for all samples, comma-separated pair " +
                "of quality threshold values for Phred and CQS quality thresholds respectively [default = 25, 30]")

def opt = cli.parse(args)

if (opt == null) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

// SEGMENTS
boolean includeNonFuncitonal = opt.'all-segments', includeAlleles = opt.'all-alleles'
if (opt.'print-library') {
    println "CDR3 extraction is possible for the following data (segments include non-functional = $includeNonFuncitonal):"
    Util.listAvailableSegments(includeNonFuncitonal, includeAlleles)
    System.exit(0)
}

if (opt.arguments().size() < 3) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

def blastPath = opt.'blast-path' ?: null

// INPUT AND OUTPUT FILES
String sampleInfoFileName = opt.'sample-metadata' ?: null, checkoutDir = opt.arguments()[0],
       assembleDir = opt.arguments()[1], outputPath = opt.arguments()[2]

if (!sampleInfoFileName && !opt.R) {
    println "[ERROR] Either --sample-metadata or --default-gene options should be provided"
    System.exit(2)
}

boolean processRaw = checkoutDir != BLANK_PATH, processAssembled = assembleDir != BLANK_PATH
if (!processRaw && !processAssembled) {
    println "[ERROR] At least one of assembled or raw (checkout) data directories should be specified"
    System.exit(2)
}
boolean processBoth = processRaw && processAssembled

// LOAD SAMPLE ENTRIES
Map rawSamplesMap = null
if (processRaw) {
    def checkoutSampleListFile = new File(checkoutDir + "/checkout.filelist.txt")
    if (!checkoutSampleListFile.exists()) {
        println "[ERROR] Sample list (checkout.filelist.txt) absent in Checkout folder $checkoutDir"
        System.exit(2)
    }

    rawSamplesMap = checkoutSampleListFile.readLines().findAll {
        !it.startsWith("#")   // no headers
    }.collectEntries {
        def splitLine = it.split('\t')
        def sampleName = splitLine[0], sampleType = splitLine[1]
        [("$sampleName\t$sampleType".toString()): splitLine[2..3]]
    }
}

Map assembledSampleMap = null
if (processAssembled) {
    def assembleSampleListFile = new File(assembleDir + "/assemble.log.txt")
    if (!assembleSampleListFile.exists()) {
        println "[ERROR] Assemble output (assemble.log.txt) absent in Assemble folder $assembleDir"
        System.exit(2)
    }

    assembledSampleMap = assembleSampleListFile.readLines().findAll {
        !it.startsWith("#")   // no headers
    }.collectEntries {
        def splitLine = it.split('\t')
        def sampleName = splitLine[0], sampleType = splitLine[1], sampleKey = "$sampleName\t$sampleType".toString()

        [(sampleKey): [splitLine[2..3], splitLine[4..5]]]
    }
}

// GENERATE SAMPLE METADATA FOR CDRBLAST
def sampleInfoMap = new HashMap()
List<String> sampleInfoLines = []

def defaultMask = opt.'default-mask' ?: "1:1", defaultSpecies = opt.'default-species' ?: "homosapiens",
    defaultFileTypes = opt.'default-file-types' ?: Util.FILE_TYPES.join(","),
    defaultQualityThreshold = opt.'default-quality-threshold' ?: "25,30", defaultChain = opt.R

if (sampleInfoFileName) {
    // Read sample metadata
    def chainsFile = new File(sampleInfoFileName)

    if (!chainsFile.exists()) {
        println "[ERROR] Sample metadata file not found"
        System.exit(2)
    }

    sampleInfoLines = chainsFile.readLines()
} else {
    // Or generate it from default parameters
    def defaultParameters = [defaultSpecies, defaultChain,
                             defaultFileTypes, defaultMask, defaultQualityThreshold].join("\t")

    (processRaw ? rawSamplesMap : assembledSampleMap).keySet().each {
        def sampleId = it.split("\t")[0]
        sampleInfoLines.add("$sampleId\t" + defaultParameters)
    }
}

// CHECK SAMPLE METADATA
sampleInfoLines.findAll { !it.startsWith("#") }.each { line ->
    def splitLine = line.split("\t")

    if (splitLine.length < 3) {
        println "[ERROR] Bad metadata line $line. Missing required fields."
        System.exit(2)
    }

    def sampleId = splitLine[0],
        species = splitLine[1], chain = splitLine[2],
        fileTypes = splitLine.length > 3 ?
                (splitLine[3] == BLANK_FIELD ? defaultFileTypes : splitLine[3]) : defaultFileTypes,
        mask = (splitLine.length > 4 ?
                (splitLine[4] == BLANK_FIELD ? defaultMask : splitLine[4]) : defaultMask),
        qualityThreshold = (splitLine.length > 5 ?
                (splitLine[5] == BLANK_FIELD ? defaultQualityThreshold : splitLine[5]) : defaultQualityThreshold).split(",")

    /*
    This will be checked at CdrBlast side
    if (!Util.isAvailable(species, chain, includeNonFuncitonal, includeAlleles)) {
        println "[ERROR] Sorry, no analysis could be performed for $species gene $chain " +
                "(include non-functional = $includeNonFuncitonal). " +
                "Possible variants are:\n"
        Util.listAvailableSegments(includeNonFuncitonal, includeAlleles)
        System.exit(-1)
    }*/
    if (!Util.MASKS.any { mask == it }) {
        println "[ERROR] Bad mask $mask on line $line. Supported masks are ${Util.MASKS.join(", ")}"
        System.exit(2)
    }
    mask.split(",").collect { it == "1" }
    fileTypes = fileTypes.split(",")
    fileTypes.each { fileType ->
        if (!Util.FILE_TYPES.any { fileType == it }) {
            println "[ERROR] Bad file type $fileType in line $line. Supported file types are ${Util.FILE_TYPES.join(", ")}"
            System.exit(2)
        }
    }

    // check if sample exists
    boolean exists = false
    fileTypes.each {
        def sampleKey = "$sampleId\t$it".toString()
        if (processRaw && rawSamplesMap.containsKey(sampleKey))
            exists = true
        if (processAssembled && assembledSampleMap.containsKey(sampleKey))
            exists = true
    }
    if (!exists) {
        println "[ERROR] Sample $sampleId with file types ${fileTypes.join(",")} does not exist in checkout/assembled sample lists"
        System.exit(2)
    }

    sampleInfoMap.put(sampleId, [species, chain, fileTypes, mask, qualityThreshold])
}

// START LOGGING
def logFile = new File(outputPath + "/cdrblast.log.txt")
if (logFile.exists())
    logFile.delete()
else
    logFile.absoluteFile.parentFile.mkdirs()

logFile.withPrintWriter { pw ->
    pw.println(Util.CDRBLAST_LOG_HEADER)

    // RUN CDRBLAST
    def baseArgs = ["--same-sample"]
    if (opt.'no-sort')
        baseArgs = [baseArgs, ["--no-sort"]]
    if (opt.p)
        baseArgs = [baseArgs, ["-p", opt.p]]
    if (blastPath)
        baseArgs = [baseArgs, ["--blast-path", blastPath]]
    if (includeNonFuncitonal)
        baseArgs = [baseArgs, ["--all-segments"]]
    String logLine
    if (processAssembled) {
        if (processBoth)
            println "[${new Date()} $scriptName] Running CdrBlast for raw and assembled data.."
        else
            println "[${new Date()} $scriptName] Running CdrBlast for assembled data.."

        sampleInfoMap.each { chainData ->
            def sampleId = chainData.key, species = chainData.value[0], chain = chainData.value[1],
                fileTypes = chainData.value[2], mask = chainData.value[3], qualityThreshold = chainData.value[4]

            def assemblyFiles = [], rawFiles = []

            fileTypes.each { fileType ->
                def sampleKey = "$sampleId\t$fileType".toString()

                if (assembledSampleMap.containsKey(sampleKey)) {
                    def correspondingFiles = assembledSampleMap[sampleKey]

                    // apply mask
                    if (fileType == 'paired') {
                        if (mask[0] + mask[1] != 0) {
                            boolean zeroMask = true
                            if (mask[0]) {
                                rawFiles.add(correspondingFiles[0][0])
                                assemblyFiles.add(correspondingFiles[1][0])
                                if (correspondingFiles[1][0] != BLANK_PATH)
                                    zeroMask = false
                            }
                            if (mask[1]) {
                                rawFiles.add(correspondingFiles[0][1])
                                assemblyFiles.add(correspondingFiles[1][1])
                                if (correspondingFiles[1][1] != BLANK_PATH)
                                    zeroMask = false
                            }
                            if (zeroMask) {
                                println "[ERROR] Sample $sampleId with mask ${mask.collect { it ? 1 : 0 }} " +
                                        "does not have any corresponding assembled files, " +
                                        "of ${correspondingFiles[1]}"
                                System.exit(2)
                            }
                        }
                    } else {
                        rawFiles.add(correspondingFiles[0][0])
                        rawFiles.add(correspondingFiles[0][1])
                        assemblyFiles.add(correspondingFiles[1][0])
                        assemblyFiles.add(correspondingFiles[1][1])
                    }
                }
            }

            // cleanup from masked files
            assemblyFiles.removeAll { it == BLANK_PATH }
            rawFiles.removeAll { it == BLANK_PATH }

            def baseArgs1 = [baseArgs, ["-R", chain], ["-S", species]]

            logLine = Util.run(new CdrBlast(), [baseArgs1, ["-a"], ["-q", qualityThreshold[1]],
                                                assemblyFiles, "$outputPath/${sampleId}.asm.cdrblast.txt"].flatten().join(" "))
            pw.println("$sampleId\t" + logLine)

            if (processBoth) {
                logLine = Util.run(new CdrBlast(), [baseArgs1, ["-q", qualityThreshold[0]],
                                                    rawFiles, "$outputPath/${sampleId}.raw.cdrblast.txt"].flatten().join(" "))
                pw.println("$sampleId\t" + logLine)
            }
        }
    } else {
        println "[${new Date()} $scriptName] Running CdrBlast for raw data.."
        sampleInfoMap.each { chainData ->
            def sampleId = chainData.key, species = chainData.value[0], chain = chainData.value[1],
                fileTypes = chainData.value[2], mask = chainData.value[3], qualityThreshold = chainData.value[4]

            def rawFiles = []

            fileTypes.each { fileType ->
                def sampleKey = "$sampleId\t$fileType".toString()

                if (rawSamplesMap.containsKey(sampleKey)) {
                    def correspondingFiles = rawSamplesMap[sampleKey]

                    // apply mask
                    if (fileType == 'paired') {
                        if (mask[0])
                            rawFiles.add(correspondingFiles[0])

                        if (mask[1])
                            rawFiles.add(correspondingFiles[1])
                    } else
                        rawFiles.add(correspondingFiles[0])
                }
            }

            logLine = Util.run(new CdrBlast(), [baseArgs, ["-R", chain], ["-S", species], ["-q", qualityThreshold[0]],
                                                rawFiles, "$outputPath/${sampleId}.raw.cdrblast.txt"].flatten().join(" "))
            pw.println("$sampleId\t" + logLine)
        }
    }
}

Util.printCmd(outputPath + "/cdrblast.cmd.txt")
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
