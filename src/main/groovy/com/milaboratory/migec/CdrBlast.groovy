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

import com.milaboratory.migec.alignment.BlastResult
import com.milaboratory.migec.clonotype.ClonotypeData
import com.milaboratory.migec.dalign.SegmentSearcher
import com.milaboratory.migec.segment.Allele
import com.milaboratory.migec.segment.Segment
import groovyx.gpars.GParsPool
import org.codehaus.groovy.runtime.ResourceGroovyMethods

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.atomic.AtomicInteger

import static com.milaboratory.migec.Util.*

//////////////
//   CLI   //
////////////
def cli = new CliBuilder(usage: "CdrBlast [options] reads1.fastq[.gz] [reads2.fastq[.gz] ...] output_file\n" +
        "NOTE: NCBI-BLAST+ package required")
cli.h("usage")
cli.R(args: 1, argName: "gene[,gene2,...]",
        "Receptor gene: 'TRA', 'TRB', 'TRG', 'TRD',  'IGL', 'IGK' or 'IGH' [required]. " +
                "Can be a comma-separated list of genes. " +
                "Use --print-library to print the list of allowed species-gene combinations.")
cli.S(args: 1, argName: "species",
        "Species: 'HomoSapiens'[default], 'MusMusculus', ... " +
                "Use --print-library to print the list of allowed species-gene combinations.")
cli.a("Input data is assembled consensuses set, not reads. " +
        "The read header must contain 'UMI:NNNNN:COUNT' entry.")
cli.q(args: 1, argName: "Phred",
        "minimum quality of read to pass filter, will be " +
                "applied only to CDR3 region [default=30 for -a or 25, out of 40]")
cli.p(args: 1,
        "number of threads to use [default = all available processors]")
cli.N(args: 1,
        "Number of reads to take. For downsampling purposes")
cli._(longOpt: "no-sort",
        "Do not sort output by clonotype count, can save some time")
cli._(longOpt: "allow-chimeric",
        "Allow V-J combinations coming from different genes")
cli._(longOpt: "blast-path", args: 1, argName: "directory",
        "Path to blast executable.")
cli._(longOpt: "same-sample", "Assembled data only (-a). Input files come from the same sample so " +
        "each UMI will be counted once.")
cli._(longOpt: "all-segments",
        "Use full V/D/J segment library (including pseudogens, etc).")
cli._(longOpt: "all-alleles",
        "Use full list of alleles (uses only major *01 alleles by default)")
cli._(longOpt: "print-library",
        "Prints out allowed species-gene pairs. " +
                "To account non-functional segment data use together with --all-segments")
cli._(longOpt: "strict-nc-handling",
        "Will not discard clonotypes with partial CDR3, " +
                "but rather try to append to corresponding complete clonotypes.")
cli._(longOpt: "debug",
        "Prints out alignment details, stores all temporary alignment files.")
cli._(longOpt: "cdr3-fastq-file", args: 1, argName: "file",
        "Store reads with extracted CDR3, appends CDR3 data in header. " +
                "Needed for 'com.milaboratory.migec.post.GroupByCdr' script.")
cli._(longOpt: "unmapped-fastq-file", args: 1, argName: "file",
        "Store reads with unmapped CDR3.")
cli._(longOpt: "cdr3-umi-table", args: 1, argName: "file",
        "Stores UMIs and their associated clonotypes. Only applicable with -a option.")
cli._(longOpt: "log-file", args: 1, argName: "file",
        "File to output cdr extraction log.")
cli._(longOpt: "log-overwrite",
        "Overwrites provided log file.")
cli._(longOpt: "log-sample-name",
        "Sample name to use in log [default = N/A].")

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

// SEGMENTS STUFF
boolean includeNonFuncitonal = opt.'all-segments', includeAlleles = opt.'all-alleles',
        allowChimeric = opt.'allow-chimeric'
if (opt.'print-library') {
    println "CDR3 extraction is possible for the following data (segments include non-functional = $includeNonFuncitonal):"
    listAvailableSegments(includeNonFuncitonal, includeAlleles)
    System.exit(0)
}

if (!opt.R) {
    println "[ERROR] Receptor gene not provided"
    cli.usage()
    System.exit(2)
}

String chainOpt = opt.R, species = opt.S ?: "HomoSapiens"

if (!chainOpt) {
    println "[ERROR] Chain argument is required for CdrBlast"
    System.exit(2)
}

def chains = chainOpt.split(",")

def hasD = chains.collect { chain ->
    if (!isAvailable(species, chain, includeNonFuncitonal, includeAlleles)) {
        println "[ERROR] Sorry, no analysis could be performed for $species gene $chain " +
                "(include non-functional = $includeNonFuncitonal). " +
                "Possible variants are:\n"
        listAvailableSegments(includeNonFuncitonal, includeAlleles)
        System.exit(2)
    }
    hasD(species, chain, includeNonFuncitonal, includeAlleles)
}

// SYSTEM
def DEBUG = opt.'debug', DEBUG_MAX_MESSAGES = 10,
    THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors(),
    SCRIPT_NAME = getClass().canonicalName

def timestamp = {
    "[${new Date()} $SCRIPT_NAME]"
}

def blastPath = opt.'blast-path' ?: ""
blastPath = blastPath.length() > 0 ? blastPath + "/" : ""

// Check for BLAST installation
try {
    ["${blastPath}convert2blastmask", "${blastPath}makeblastdb", "${blastPath}blastn"].each { it.execute().waitFor() }
} catch (IOException e) {
    println "[ERROR] Problems with BLAST installation. " + e.message
    System.exit(2)
}

// INPUT, OUTPUT AND TEMPORARY FILES
if (opt.arguments().size() < 2) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

def inputFileNames = opt.arguments()[0..-2].collect { it.toString() },
    outputFileName = opt.arguments()[-1].toString()

def TMP_FOLDER = new File(inputFileNames[0]).absolutePath + "-cdrblast-" + UUID.randomUUID().toString()

def TMP_FOLDER_FILE = new File(TMP_FOLDER)
TMP_FOLDER_FILE.mkdirs()
//if (!DEBUG)
//    TMP_FOLDER_FILE.deleteOnExit()

def outputFile = null
if (outputFileName != BLANK_PATH && (outputFile = new File(outputFileName)).parentFile)
    outputFile.parentFile.mkdirs()

// BLAST SETTINGS
int ALLELE_SEED_INNER = 10, ALLELE_SEED_OUTER = 6,
    ALLELE_TAIL_V_MAX = 50, ALLELE_TAIL_J_MAX = 20,
    GAP_OPEN = 5, GAP_EXTEND = 2, WORD_SIZE = 7, REWARD = 2, PENALTY = -3

String BLAST_FLAGS = "-lcase_masking"

// D ALIGNER SETTINGS
int MIN_D_MATCH = 2

// BLAST results filtering
int TOP_SEQS = 1,
    MIN_SEGMENT_IDENT = 7,
    MIN_CDR_LEN = 7, // at least 1 nt + conserved AAs
    MAX_CDR_LEN = 300 // there are some IGH with 35+ AAs in CDR3

// LOGGING
String logFileName = opt.'log-file' ?: null
boolean overwriteLog = opt.'log-overwrite'
String sampleName = opt.'log-sample-name' ?: "N/A"

// PARAMETERS
def cdr3FastqFile = opt.'cdr3-fastq-file', cdr3UmiTable = opt.'cdr3-umi-table',
    unmappedFastqFile = opt.'unmapped-fastq-file'

boolean strictNcHandling = opt.'strict-nc-handling',
        assembledInput = opt.a,
        doSort = !opt.'no-sort',
        sameSample = opt.'same-sample'

int qualThreshold = opt.q ? Integer.parseInt(opt.q) : (opt.a ? 30 : 25),
    nReads = Integer.parseInt(opt.N ?: "-1")

/////////////////////////////////////////////////////////////
// Stage 0: Pre-load segment information. Set up aligners //
///////////////////////////////////////////////////////////

def segments = new HashMap<String, Segment>(), alleles = new HashMap<String, Allele>()
def vAlleles = new ArrayList<Allele>(), jAlleles = new ArrayList<Allele>(),
    dAlleles = new ArrayList<Allele>()
def collapseAlleleMap = new HashMap<String, Allele>()

chains.eachWithIndex { String chain, int ind ->
    println "${timestamp()} Loading ${chain} Variable and Joining segment data"

    def resFile = getSegmentsFile(includeNonFuncitonal, includeAlleles)

    resFile.splitEachLine("\t") {
        if (species.toUpperCase() == it[0].toUpperCase() &&
                chain.toUpperCase() == it[1].toUpperCase()) { // Take only alleles of a given chain and species
            def type = it[2].charAt(0).toUpperCase(), seq = it[5], refPoint = Integer.parseInt(it[4])

            if (type == "D") {
                // Simply store D alleles
                dAlleles.add(new Allele(alleleId: it[3],
                        segmentId: it[3].split("[*]")[0],
                        seq: seq))
                // note that D segment could have different strand
                dAlleles.add(new Allele(alleleId: it[3],
                        segmentId: it[3].split("[*]")[0],
                        seq: revCompl(seq)))
            } else {
                // More complex processing for V/J alleles
                // Mask all except seed region
                seq = seq.toLowerCase()
                seq = seq.toCharArray()
                int tailFrom = (type == "V" ? ALLELE_SEED_OUTER : ALLELE_SEED_INNER),
                    tailTo = (type == "V" ? ALLELE_SEED_INNER : ALLELE_SEED_OUTER)
                for (int i = -tailFrom; i < tailTo; i++)
                    if (refPoint + i >= 0 && refPoint + i < seq.length)
                        seq[refPoint + i] = seq[refPoint + i].toUpperCase()
                seq = new String(seq)

                // Trim extending tails of alleles which are not covered by read
                if (type == "V") {
                    if (ALLELE_TAIL_V_MAX <= refPoint) {
                        seq = seq.substring(refPoint - ALLELE_TAIL_V_MAX)
                        refPoint = ALLELE_TAIL_V_MAX
                    }
                } else
                    seq = seq.substring(0, Math.min(refPoint + ALLELE_TAIL_J_MAX + 1, seq.length()))

                def allele = collapseAlleleMap[seq] // Collapse alleles with identical sequences
                boolean newAllele = allele == null
                if (newAllele) {
                    // it[0] - species
                    allele = new Allele(alleleId: it[3], segmentId: it[3].split("[*]")[0],
                            refPoint: refPoint,
                            seq: seq,
                            type: type,
                            chain: chain,
                            hasD: hasD[ind])
                    collapseAlleleMap.put(seq, allele)
                } else
                    allele.alleleId += "," + it[3]

                if (newAllele) {
                    // V and J alleles separately - for alignment code simplicity
                    if (allele.type == "V")
                        vAlleles.add(allele)
                    else
                        jAlleles.add(allele)

                    // Segment info - for allele frequencies
                    def segment = segments[allele.segmentId]
                    if (segment == null)
                        segment = new Segment(id: allele.segmentId, type: allele.type)
                    segment.alleles.add(allele)
                    segments.put(segment.id, segment)
                }
            }
        }
    }

    collapseAlleleMap.values().each {
        alleles.put(it.alleleId, it)
    }
}

// Set up segment alignment
// Make blast db
def chainSign = chains.join("+")
["V", "J"].each { seg ->
    println "${timestamp()} Generating BLAST database for $chainSign $seg"
    new File("$TMP_FOLDER/${chainSign}_${seg}.fa").withPrintWriter { pw ->
        alleles.each { allele ->
            if (allele.value.type == seg)
                pw.println(">$allele.key\n$allele.value.seq")
        }
    }

    ("${blastPath}convert2blastmask -in $TMP_FOLDER/${chainSign}_${seg}.fa -out $TMP_FOLDER/${chainSign}_${seg}.msk " +
            "-masking_algorithm 'CDRBLAST' -masking_options 'NA'").execute().waitFor()

    ("${blastPath}makeblastdb -in $TMP_FOLDER/${chainSign}_${seg}.fa -mask_data $TMP_FOLDER/${chainSign}_${seg}.msk " +
            "-dbtype nucl -out $TMP_FOLDER/${chainSign}_${seg}").execute().waitFor()
}

// D segment searcher
def dSegmentSearcher = new SegmentSearcher(dAlleles, MIN_D_MATCH)

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Stage I: blast for unique sequences. Here we ONLY determine V and J segment presence and their IDs //
// seqSet (unique sequences) -> blast mapping data (unique sequence + V and J mapping)               //
//////////////////////////////////////////////////////////////////////////////////////////////////////
def seqMap = new HashMap<String, Integer>(),
    seqList = new ArrayList<String>()
int counter = 0, uniqueCounter = 0
inputFileNames.each { inputFileName ->
    println "${timestamp()} Reading $inputFileName and generating list of unique sequences to map V and J genes"
    def reader = getReader(inputFileName)
    String header
    while (((header = reader.readLine()) != null) && (nReads < 0 || counter < nReads)) {
        if (assembledInput) {
            if (!header.contains("UMI:")) {
                println "[ERROR] Assembled input specified, but no UMI field in header"
                System.exit(2)
            } else if (!header.split("[@ ]").find { it.startsWith("UMI") }.split(":")[2].isInteger()) {
                println "[ERROR] Assembled input specified, but UMI field in header does not contain count, " +
                        "maybe you're running on raw checkout output?"
                System.exit(2)
            }
        }

        def seq = reader.readLine()
        if (!seqMap.containsKey(seq)) {
            seqMap.put(seq, uniqueCounter) // build a set of unique sequences to reduce the amount of data
            seqList.add(seq)
            uniqueCounter++
        }
        counter++
        reader.readLine()
        reader.readLine()

        if (counter % 1000000 == 0)
            println "${timestamp()} $counter reads processed, $uniqueCounter unique"
    }
}

def queryFilePrefix = TMP_FOLDER + "/" + getFastqPrefix(inputFileNames[0])
println "[${new Date()} $SCRIPT_NAME] Making a temporary FASTA file for BLAST queries"
THREADS = Math.min(THREADS, seqList.size())
int chunk_sz = seqList.size() / THREADS
for (int p = 0; p < THREADS; p++) { // split fasta for blast parallelizaiton
    def file = new File(queryFilePrefix + "_${p}.fa")
    int end = p == THREADS - 1 ? seqList.size() : (p + 1) * chunk_sz
    end = Math.min(end, seqList.size())
    file.withPrintWriter { pw ->
        for (int i = p * chunk_sz; i < end; i++) {
            pw.println(">$i\n${seqList[i]}")
        }
    }
    //if (!DEBUG)
    //    file.deleteOnExit()
}

["V", "J"].each { seg -> // run blast for v and j segments separately
    println "${timestamp()} Pre-aligning $chainSign $seg segment with BLAST <$THREADS threads>"
    def blastProcs = (0..(THREADS - 1)).collect { p ->
        def blastOutFname = "${queryFilePrefix}_${chainSign}_${seg}_${p}.blast" // temp
        //if (!DEBUG)
        //    new File(blastOutFname).deleteOnExit()

        // A trick to pass -outfmt argument correctly
        def blastCmd = ["${blastPath}blastn",
                        "-query", "${queryFilePrefix}_${p}.fa",
                        BLAST_FLAGS.split(" "),
                        "-gapopen", "$GAP_OPEN", "-gapextend", "$GAP_EXTEND",
                        "-word_size", "$WORD_SIZE", "-reward", "$REWARD", "-penalty", "$PENALTY",
                        "-db", "$TMP_FOLDER/${chainSign}_${seg}",
                        "-max_target_seqs", "$TOP_SEQS", "-out", "$blastOutFname"]

        blastCmd = [blastCmd, "-outfmt", "6 qseqid sseqid qstart qend sstart send qseq sseq nident"].flatten()
        blastCmd.execute()
        //"bash $TMP_FOLDER/runblast.sh ${queryFilePrefix}_${p}.fa $TMP_FOLDER/${chainSign}_${seg} $blastOutFname".execute()
    }

    blastProcs.each { it.waitFor() }  // This is the only way to parallelize blast
    blastProcs.each { it.closeStreams() }
    blastProcs.each { it.destroy() }

    println "[${new Date()} $SCRIPT_NAME] Done"
}

//////////////////////////////////////////////////
// Stage II: Read blast results & extract CDRs //
////////////////////////////////////////////////

def vMappingData = new HashMap<Integer, BlastResult>(),
    jMappingData = new HashMap<Integer, BlastResult>()

println "${timestamp()} Reading in BLAST results"
["V", "J"].eachWithIndex { segment, ind -> // parse blast output
    (0..<THREADS).each { p ->
        new File("${queryFilePrefix}_${chainSign}_${segment}_${p}.blast").splitEachLine("\t") { line ->
            int nIdent = Integer.parseInt(line[8])
            if (nIdent >= MIN_SEGMENT_IDENT) {
                def id = Integer.parseInt(line[0])
                def result = new BlastResult(
                        alleles[line[1]],
                        Integer.parseInt(line[2]) - 1, Integer.parseInt(line[3]) - 1,
                        Integer.parseInt(line[4]) - 1, Integer.parseInt(line[5]) - 1,
                        line[6],
                        nIdent
                )
                def mappingData = segment == 'V' ? vMappingData : jMappingData
                def prevMapping = mappingData[id]

                if (prevMapping == null || prevMapping.nIdent < nIdent)
                    mappingData.put(id, result)   // take result with highest identity
            }
        }
    }
}
println "${timestamp()} Done ${vMappingData.size()} V and ${jMappingData.size()} J mapped, " +
        "out of total ${seqMap.size()} sequence variants"


def readId2ClonotypeData = new ConcurrentHashMap<Integer, ClonotypeData>()
def uniqueClonotypes = new ConcurrentHashMap<ClonotypeData, ClonotypeData>() // util
println "${timestamp()} Extracting CDRs"

def vRefNotInAlign = new AtomicInteger(0), jRefNotInAlign = new AtomicInteger(0),
    reverseFound = new AtomicInteger(0), badOrientation = new AtomicInteger(0),
    cdrStartOut = new AtomicInteger(0), cdrEndOut = new AtomicInteger(0),
    shortCdr = new AtomicInteger(0), longCdr = new AtomicInteger(0)
def badExamples = new AtomicInteger(0), goodExamples = new AtomicInteger(0)

def debugMessagesGood = Collections.synchronizedList(new ArrayList<String>()),
    debugMessagesBad = Collections.synchronizedList(new ArrayList<String>())

GParsPool.withPool THREADS, {
    seqMap.eachParallel { entry ->
        String seq = entry.key
        int id = entry.value
        def vMapping = vMappingData[id], jMapping = jMappingData[id]
        if (vMapping != null && jMapping != null &&
                (allowChimeric || vMapping.allele.chain == jMapping.allele.chain)) {
            Allele vAllele = vMapping.allele,
                   jAllele = jMapping.allele
            boolean vRC = vMapping.aFrom > vMapping.aTo, jRC = jMapping.aFrom > jMapping.aTo

            // In case of rc flip so that start<end
            if (vRC) {
                int temp = vMapping.aTo
                vMapping.aTo = vMapping.aFrom
                vMapping.aFrom = temp
            }

            if (jRC) {
                int temp = jMapping.aTo
                jMapping.aTo = jMapping.aFrom
                jMapping.aFrom = temp
            }

            int vRef = vAllele.refPoint,
                jRef = jAllele.refPoint

            boolean vRefInAlign = vRef + 1 >= vMapping.aFrom && vRef - 1 <= vMapping.aTo,
                    jRefInAlign = jRef + 1 >= jMapping.aFrom && jRef - 1 <= jMapping.aTo

            boolean bad = true

            int cdrFrom, cdrTo, cdrLen

            if (vRefInAlign && jRefInAlign && vRC == jRC) {
                if (vRC) { // after that step RC is same to non-RC
                    reverseFound.get()
                    seq = revCompl(seq)

                    int temp = vMapping.qTo
                    vMapping.qTo = seq.length() - vMapping.qFrom - 1
                    vMapping.qFrom = seq.length() - temp - 1

                    temp = jMapping.qTo
                    jMapping.qTo = seq.length() - jMapping.qFrom - 1
                    jMapping.qFrom = seq.length() - temp - 1
                }

                // Determine v&j ref positions with respect to gaps
                vRef -= vMapping.aFrom
                int fixedRef = vRef
                for (int i = 0; i < fixedRef; i++)
                    if (vMapping.gaps[i])
                        vRef--
                vRef += vMapping.qFrom

                jRef -= jMapping.aFrom
                fixedRef = jRef
                for (int i = 0; i < fixedRef; i++)
                    if (jMapping.gaps[i])
                        jRef--
                jRef += jMapping.qFrom

                // Extract CDR3
                // Cys & Phe are included
                cdrFrom = vRef - 3
                cdrTo = jRef + 4
                cdrLen = cdrTo - cdrFrom

                if (cdrFrom >= 0 && cdrTo <= seq.length() && cdrLen >= MIN_CDR_LEN && cdrLen <= MAX_CDR_LEN) {
                    // Do not create new CDR objects - save memory
                    // CDR object has to store sequence string for hashing
                    def clonotypeData = new ClonotypeData(vAllele, jAllele,
                            seq.substring(cdrFrom, cdrTo),
                            cdrFrom, cdrTo, vMapping.qTo, jMapping.qFrom, vRC) // to get hash

                    clonotypeData =
                            uniqueClonotypes.putIfAbsent(clonotypeData, clonotypeData) ?: // use already generated
                                    clonotypeData // create new

                    // Unique read identifier -> CDR3 extraction result
                    readId2ClonotypeData.put(id, clonotypeData)

                    bad = false
                } else {
                    if (cdrFrom < 0)
                        cdrStartOut.incrementAndGet()
                    if (cdrTo > seq.length())
                        cdrEndOut.incrementAndGet()
                    if (cdrLen < MIN_CDR_LEN)
                        shortCdr.incrementAndGet()
                    if (cdrLen > MAX_CDR_LEN)
                        longCdr.incrementAndGet()
                }
            } else {
                if (!vRefInAlign)
                    vRefNotInAlign.incrementAndGet()
                if (!jRefInAlign)
                    jRefNotInAlign.incrementAndGet()
                if (vRC != jRC)
                    badOrientation.incrementAndGet()
            }

            if (DEBUG) {
                if (bad && badExamples.get() < DEBUG_MAX_MESSAGES) {
                    def message = ((vRC == jRC && vRC) ? revCompl(seq) : seq) + "\n"

                    if (!vRefInAlign)
                        message += vAllele.alleleId + "\t" +
                                (vRC ? revCompl(vAllele.seq.toUpperCase()) : vAllele.seq.toUpperCase()) + "\t" +
                                "aFrom=$vMapping.aFrom aTo=$vMapping.aTo ref=$vRef\n"

                    if (!jRefInAlign)
                        message += jAllele.alleleId + "\t" +
                                (vRC ? revCompl(jAllele.seq.toUpperCase()) : jAllele.seq.toUpperCase()) + "\t" +
                                "aFrom=$jMapping.aFrom aTo=$jMapping.aTo ref=$jRef\n"

                    debugMessagesBad.add(message)

                    badExamples.incrementAndGet()
                } else if (goodExamples.get() < DEBUG_MAX_MESSAGES) {
                    def message = seq + "\n"// already reversed

                    message += seq.substring(cdrFrom, cdrTo) + "\n"

                    message += vAllele.alleleId + "\t" +
                            vAllele.seq.toUpperCase() + "\t" +
                            "qTo=$vMapping.qTo qFrom=$vMapping.qFrom aTo=$vMapping.aTo aFrom=$vMapping.aFrom ref=$vRef\n"

                    message += jAllele.alleleId + "\t" + revCompl(jAllele.seq.toUpperCase()) + "\t" +
                            "qTo=$jMapping.qTo qFrom=$jMapping.qFrom aTo=$jMapping.aTo aFrom=$jMapping.aFrom ref=$jRef\n"

                    debugMessagesGood.add(message)

                    goodExamples.incrementAndGet()
                }
            }
        }
    }
}

if (DEBUG) {
    println "Reporting $DEBUG_MAX_MESSAGES good and bad alignments"
    debugMessagesGood.eachWithIndex { it, ind -> print ">GOOD$ind\n$it" }
    debugMessagesBad.eachWithIndex { it, ind -> print ">BAD$ind\n$it" }
}

println "${timestamp()} Done. ${readId2ClonotypeData.size()} CDRs mapped, " +
        "out of total ${seqMap.size()} sequence variants"

if (DEBUG) {
    println "reverseFound=" + reverseFound.get() + " (${(int) (reverseFound.get() / seqMap.size() * 100)}%)"
    println "badOrientation=" + badOrientation.get()
    println "vRefNotInAlignment=" + vRefNotInAlign.get()
    println "jRefNotInAlignment=" + jRefNotInAlign.get()
    println "cdrStartOut=" + cdrStartOut.get()
    println "cdrEndOut=" + cdrEndOut.get()
    println "shortCdr3=" + shortCdr.get() + " (<$MIN_CDR_LEN nt)"
    println "longCdr3=" + longCdr.get() + " (>$MAX_CDR_LEN nt)"
}

// Mapping Diversity segments
if (!dAlleles.empty) {
    println "${timestamp()} Done. Mapping Diversity segments"
    GParsPool.withPool THREADS, {
        uniqueClonotypes.values().eachParallel { ClonotypeData clonotypeData ->
            clonotypeData.appendDInfo(dSegmentSearcher)
        }
    }
}

//////////////////////////////////
// Stage III: Append read data //
////////////////////////////////
counter = 0
int goodReads = 0, goodEvents = 0, mappedReads = 0, mappedEvents = 0, totalReads = 0, totalEvents = 0
def writerCdr3Fastq = cdr3FastqFile ? getWriter(cdr3FastqFile) : null,
    writerUnmapped3Fastq = unmappedFastqFile ? getWriter(unmappedFastqFile) : null,
    writerCdr3Umi = assembledInput && cdr3UmiTable ? getWriter(cdr3UmiTable) : null

if (writerCdr3Umi)
    writerCdr3Umi.println("#umi\tmig_sz\tcdr3nt\tv\tj\td")

def usedUmis = new HashSet<String>()

inputFileNames.each { inputFileName ->
    println "${timestamp()} Appending read data from $inputFileName"
    def reader = getReader(inputFileName)
    while (((header = reader.readLine()) != null) && (nReads < 0 || counter < nReads)) {
        def seq = reader.readLine()
        reader.readLine()
        def qual = reader.readLine(), oldQual = qual
        int id = seqMap[seq]  // get unique read sequence identifier
        def clonotypeData = readId2ClonotypeData[id] // corresponding CDR3 extraction result

        int increment = 1
        String umi
        boolean duplicateUmi = false

        if (assembledInput) {
            def splitHeader = header.split("[@ ]")
            def umiFields = splitHeader.find { it.startsWith("UMI") }.split(":")
            umi = umiFields[1]
            increment = Integer.parseInt(umiFields[2])
        }

        // Has CDR3
        if (clonotypeData != null) {
            if (clonotypeData.rc) // found in RC
                qual = qual.reverse()

            qual = qual.substring(clonotypeData.cdrFrom, clonotypeData.cdrTo)

            // Increment counters if quality is good
            if (minQual(qual) >= qualThreshold) {
                if (assembledInput && sameSample) {
                    def umiCdr = umi + "_" + clonotypeData.cdr3Seq
                    duplicateUmi = usedUmis.contains(umi)
                    if (!duplicateUmi)
                        usedUmis.add(umiCdr)
                }

                if (!duplicateUmi) {
                    // Protect from duplicate counting when UMI overlap
                    clonotypeData.counts[0]++
                }
                goodEvents++

                clonotypeData.counts[2] += increment
                goodReads += increment

                if (cdr3FastqFile)
                    writerCdr3Fastq.writeLine(header + " CDR3:" + clonotypeData.cdr3Seq + ":" +
                            clonotypeData.cdrFrom +
                            " V:" + clonotypeData.vAllele.alleleId +
                            " J:" + clonotypeData.jAllele.alleleId +
                            " D:" + clonotypeData.dAllele +
                            "\n" +
                            seq + "\n+\n" + oldQual)

                if (writerCdr3Umi)
                    writerCdr3Umi.println([umi, increment,
                                           clonotypeData.cdr3Seq,
                                           clonotypeData.vAllele.alleleId,
                                           clonotypeData.jAllele.alleleId,
                                           clonotypeData.dAllele].join("\t"))
            }

            // Total (good+bad) for cases in which CDR3 was extracted
            if (!duplicateUmi) {
                clonotypeData.counts[1]++
            }
            mappedEvents++
            clonotypeData.counts[3] += increment
            mappedReads += increment
        } else if (unmappedFastqFile) {
            writerUnmapped3Fastq.println(header + "\n" + seq + "\n+\n" + qual)
        }

        // Total, including cases when CDR3 was not extracted
        totalEvents++
        totalReads += increment
        counter++

        if (counter % 1000000 == 0)
            println "[${new Date()} $SCRIPT_NAME] $counter reads processed"
    }
}
if (cdr3FastqFile)
    writerCdr3Fastq.close()
if (writerCdr3Umi)
    writerCdr3Umi.close()
if (unmappedFastqFile)
    writerUnmapped3Fastq.close()

println "${timestamp()} Done"
println "EVENTS (good mapped total):\t$goodEvents\t$mappedEvents\t$totalEvents\t|\t" +
        "${(int) (100 * (double) goodEvents / (double) totalEvents)}%\t" +
        "${(int) (100 * (double) mappedEvents / (double) totalEvents)}%\t100%"
println "READS (good mapped total):\t$goodReads\t$mappedReads\t$totalReads\t|\t" +
        "${(int) (100 * (double) goodReads / (double) totalReads)}%\t" +
        "${(int) (100 * (double) mappedReads / (double) totalReads)}%\t100%"

// Combine clonotypes based on signature, not shortSignature
// shortSignature used earlier needed for backtracking reads and is redundant because of cdrFrom
def cloneMap = new HashMap<String, int[]>()
def nonCanonicalCounts = new HashMap<String, int[]>()

if (strictNcHandling)
    uniqueClonotypes.keySet().each {
        if (!it.isCanonical()) {
            def counts = nonCanonicalCounts.get(it.cdr3Seq)
            if (counts == null)
                nonCanonicalCounts.put(it.cdr3Seq, counts = new int[4])
            for (int i = 0; i < 4; i++)
                counts[i] += it.counts[i]
        }
    }

uniqueClonotypes.keySet().each {
    if (it.isCanonical()) {
        def signature = it.signature
        def counts = cloneMap.get(signature)
        if (counts == null)
            cloneMap.put(signature, counts = new int[4])
        for (int i = 0; i < 4; i++)
            counts[i] += it.counts[i]

        if (strictNcHandling) {
            def ncToRemove = new HashSet<String>()
            nonCanonicalCounts.entrySet().each { ncEntry ->
                if (ncEntry.key.contains(it.cdr3Seq)) {
                    for (int i = 0; i < 4; i++)
                        counts[i] += ncEntry.value[i]
                    ncToRemove.add(ncEntry.key)
                }
            }
            ncToRemove.each { nc2remove ->
                nonCanonicalCounts.remove(nc2remove)
            }
        }
    }
}

if (outputFile) {
    println "${timestamp()} Writing output"
    outputFile.withPrintWriter { pw ->
        pw.println("Count\tFraction\t" +
                "CDR3 nucleotide sequence\tCDR3 amino acid sequence\t" +
                "V segments\tJ segments\tD segments\t" +
                "Last V nucleotide position\t" +
                "First D nucleotide position\tLast D nucleotide position\t" +
                "First J nucleotide position\t" +
                "Good events\tTotal events\tGood reads\tTotal reads")

        (doSort ? cloneMap.sort { -it.value[0] } : cloneMap).each {
            if (it.value[0] > 0)
                pw.println(it.value[0] + "\t" + (it.value[0] / goodEvents) + // percentage from good events only
                        "\t" + it.key +
                        "\t" + it.value.collect().join("\t"))
        }
    }
}

// Those were created by blast and have to be removed manually
if (!DEBUG)
    ResourceGroovyMethods.deleteDir(TMP_FOLDER_FILE)
//TMP_FOLDER_FILE.listFiles().each { it.deleteOnExit() }

// Append to log and report to batch runner
def logLine = [(assembledInput ? "asm" : "raw"), outputFile ? outputFile.absolutePath : BLANK_PATH,
               inputFileNames.join(","),
               goodEvents, mappedEvents, totalEvents,
               goodReads, mappedReads, totalReads].join("\t")

if (logFileName) {
    def logFile = new File(logFileName)

    if (logFile.exists()) {
        if (overwriteLog)
            logFile.delete()
    } else {
        logFile.absoluteFile.parentFile.mkdirs()
        logFile.withPrintWriter { pw ->
            pw.println(CDRBLAST_LOG_HEADER)
        }
    }

    logFile.withWriterAppend { logWriter ->
        logWriter.println("$sampleName\t" + logLine)
    }
}

return logLine
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
