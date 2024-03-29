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

import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.atomic.AtomicLong
import java.util.regex.Pattern

import static com.milaboratory.migec.Util.BLANK_PATH

def mm = "15:0.2:0.05", rcm = "0:1", mtrim = "10", omf = "5", oss = "5", oms = "10"
def cli = new CliBuilder(usage:
        "Checkout [options] barcode_file R1.fastq[.gz] [R2.fastq[.gz] or ${BLANK_PATH}] output_dir/'")
cli.h("usage")
cli.o('Oriented reads, so master barcode has to be in R1. ' +
        'Default: scans both reads (if pair-end) for master barcode')
cli.u('Save UMI region specified by capital N\'s in barcode sequence to the header')
cli.t('Trim barcode sequences. Will also trim boundary sequence between barcode and read start/end if it is short enough.')
cli.e('Remove template-switching trace, ^T{0,3}G{3,7}. Only works in compination with -t')
cli._(longOpt: 'max-trim-nts', args: 1, argName: "int",
        'Maximal number of nucleotides between barcode and read start/end that could be trimmed')
cli.r(args: 1, "RC mask to apply after master-slave is determined, e.g. will RC slave read if set to 0:1. " +
        "Default: $rcm")
cli.p(args: 1, 'Number of threads. Default: all available processors.')
cli._(args: 1, longOpt: 'first', 'Number of reads to try. If <0 will process all reads. [default = -1]')
cli._(longOpt: 'overlap', 'Will try to overlap paired reads.')
cli._(longOpt: 'overlap-max-offset', args: 1, "Max offset for overlap. " +
        "Shoud be increased if reads are read-through. [default = $omf]")
cli._(longOpt: 'overlap-seed-size', args: 1,
        "Number of perfectly matching nucleotides to perform overlap, should be. [default = $oss]")
cli._(longOpt: 'overlap-match-size', args: 1,
        "Number of fuzzy matching nucleotides to perform overlap, should be. [default = $oms]")
cli._(longOpt: 'overlap-allow-partial',
        "Allow partial fuzzy overlap due to reaching end of one of reads.")
cli._(longOpt: 'rc-barcodes', 'Also search for reverse-complement barcodes.')
cli._(longOpt: 'skip-undef', 'Do not store reads that have no barcode match in separate file.')
cli.m(args: 1,
        argName: 'LQ:E0:E1', "Low quality threshold : Ratio of low-quality errors : Ratio of high-quality errors. " +
        "Used to calculate the maximum number of allowed mismatches when performing alignment to full barcode, " +
        "ratios are scaled to full barcode length. Default: $mm")
cli._(longOpt: 'append', 'Will append to existing directory. ' +
        'Useful when samples were previously split using Illumina indices')
cli.c('Compressed output')
def opt = cli.parse(args)
if (opt == null || opt.arguments().size() < 4) {
    println "[ERROR] Too few arguments provided"
    cli.usage()
    System.exit(2)
}

if (opt.h) {
    cli.usage()
    System.exit(0)
}

def scriptName = getClass().canonicalName
boolean compressed = opt.c, oriented = opt.o, extractUMI = opt.u, addRevComplBc = opt."rc-barcodes",
        trimBc = opt.t, removeTS = opt.e, overlap = opt.'overlap', noUndef = opt.'skip-undef',
        append = opt.'append'
int trimSizeThreshold = (opt.'max-trim-nts' ?: mtrim).toInteger()
int firstReadsToTake = (opt.'first' ?: "-1").toInteger()
int maxOverlapOffset = (opt.'overlap-max-offset' ?: omf).toInteger(),
    overlapFuzzySize = (opt.'overlap-match-size' ?: oms).toInteger(),
    overlapSeedSize = (opt.'overlap-seed-size' ?: oss).toInteger()
boolean allowPartialOverlap = opt.'overlap-allow-partial'
def rcReadMask = (opt.r ?: rcm).split(":").collect { Integer.parseInt(it) > 0 }
int THREADS = opt.p ? Integer.parseInt(opt.p) : Runtime.getRuntime().availableProcessors()
def barcodesFileName = opt.arguments()[0],
    inputFileName1 = opt.arguments()[1],
    inputFileName2 = opt.arguments()[2],
    outputDir = opt.arguments()[3]
def mmData = (opt.m ?: mm).split(":").collect { Double.parseDouble(it) }
def paired = inputFileName2 != BLANK_PATH

if (overlap && !paired) {
    println "[ERROR] Overlap requested for unpaired reads"
    System.exit(2)
}

//========================
//  MISC IUPAC TABLES
//========================
def complements = ['A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'y': 'r',
                   'r': 'y', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n']

def redundancy = ['A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': '[ATGC]', 'R': '[AG]', 'Y': '[CT]', 'M': '[AC]',
                  'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'V': '[ACG]', 'D': '[AGT]', 'H': '[ACT]', 'B': '[CGT]']

def redundancy2 = new HashMap<Character, Set<Character>>([
        ((Character) 'A'): new HashSet<Character>([(Character) 'A']),
        ((Character) 'T'): new HashSet<Character>([(Character) 'T']),
        ((Character) 'G'): new HashSet<Character>([(Character) 'G']),
        ((Character) 'C'): new HashSet<Character>([(Character) 'C']),
        ((Character) 'N'): new HashSet<Character>([(Character) 'A', (Character) 'T', (Character) 'G', (Character) 'C']),
        ((Character) 'R'): new HashSet<Character>([(Character) 'A', (Character) 'G']),
        ((Character) 'Y'): new HashSet<Character>([(Character) 'C', (Character) 'T']),
        ((Character) 'M'): new HashSet<Character>([(Character) 'A', (Character) 'C']),
        ((Character) 'S'): new HashSet<Character>([(Character) 'G', (Character) 'C']),
        ((Character) 'W'): new HashSet<Character>([(Character) 'A', (Character) 'T']),
        ((Character) 'K'): new HashSet<Character>([(Character) 'G', (Character) 'T']),
        ((Character) 'V'): new HashSet<Character>([(Character) 'A', (Character) 'C', (Character) 'G']),
        ((Character) 'D'): new HashSet<Character>([(Character) 'A', (Character) 'G', (Character) 'T']),
        ((Character) 'H'): new HashSet<Character>([(Character) 'A', (Character) 'C', (Character) 'T']),
        ((Character) 'B'): new HashSet<Character>([(Character) 'C', (Character) 'G', (Character) 'T'])])

//========================
//       BARCODE DATA
//========================
def sampleIds = new ArrayList<String>()
def barcodeList = new List<String>[2]
def seeds = new List<Pattern>[2]
def umiPositions = new List<List<Integer>>[2]
def maxGoodMMs = new List<Integer>[2], maxBadMMs = new List<Integer>[2]

// Init
barcodeList[0] = new ArrayList<String>()
barcodeList[1] = new ArrayList<String>()
seeds[0] = new ArrayList<Pattern>()
seeds[1] = new ArrayList<Pattern>()
umiPositions[0] = new ArrayList<List<Integer>>()
umiPositions[1] = new ArrayList<List<Integer>>()
maxBadMMs[0] = new ArrayList<Integer>()
maxBadMMs[1] = new ArrayList<Integer>()
maxGoodMMs[0] = new ArrayList<Integer>()
maxGoodMMs[1] = new ArrayList<Integer>()

// Subs for parsing barcode data
def getSeed = { String barcode ->
    Pattern.compile(barcode.replaceAll(/./) { redundancy."$it" ?: '[ATGC]' }) // redundancy has only uppecase
}

def parseUmiPositions = { String barcode ->
    def posList = new ArrayList<Integer>()
    for (int i = 0; i < barcode.length(); i++)
        if (barcode.charAt(i) == 'N')
            posList.add(i)
    posList
}

def addBarcode = { String barcode, int slave ->
    barcodeList[slave].add(barcode.toUpperCase()) // important - no lcase masking later
    seeds[slave].add(getSeed(barcode))
    umiPositions[slave].add(parseUmiPositions(barcode))
    maxBadMMs[slave].add((int) (barcode.length() * mmData[1]))
    maxGoodMMs[slave].add((int) (barcode.length() * mmData[2]))
}

// Load barcode data from file
new File(barcodesFileName).splitEachLine("[\t ]") { sl ->
    if (sl.size() > 0) {
        sl = sl.collect { it.length() == 0 ? BLANK_PATH : it }
        // filter by input file name
        if (
        (inputFileName2 == BLANK_PATH && (sl.size() < 4 || sl[3] == BLANK_PATH || (inputFileName1 == sl[3]))) ||
                (sl.size() < 5 ||
                        ((sl[3] == BLANK_PATH || inputFileName1 == sl[3]) && (sl[4] == BLANK_PATH || inputFileName2 == sl[4])))
        ) {
            if (!sl[0].startsWith("#")) {
                sampleIds.add(sl[0])
                if (addRevComplBc)
                    sampleIds.add(sl[0])
                if (sl[1] != null && (sl[1] = sl[1].trim()).length() > 0 &&
                        sl[1].toCharArray().every { complements.keySet().contains((String) it) }) {
                    addBarcode(sl[1], 0)
                    if (addRevComplBc)
                        addBarcode(Util.revComplExt(sl[1]), 0)
                } else {
                    println "[ERROR] Bad barcode ${sl[1]}. Terminating"
                    System.exit(2)
                }
                if (sl[2] != null && (sl[2] = sl[2].trim()).length() > 0 &&
                        sl[2].toCharArray().every { complements.keySet().contains((String) it) }) {
                    addBarcode(sl[2], 1)
                    if (addRevComplBc)
                        addBarcode(Util.revComplExt(sl[2]), 1)
                } else {
                    addBarcode('n', 1)
                    if (addRevComplBc)
                        addBarcode('n', 1)
                }
            }
        }
    }
}

println "[${new Date()} $scriptName] Loaded barcode data for " +
        "${addRevComplBc ? (sampleIds.size() / 2) : sampleIds.size()} sample(s)"

//========================
//  BARCODE SEARCH UTILS
//========================
def findAllMatches = { String seq, Pattern seed ->  // all matches for seed
    def positions = new ArrayList<Integer>()
    def matcher = seed.matcher(seq)
    while (matcher.find()) {
        positions.add(matcher.start())
    }
    positions
}
def lq = mmData[0]
def hasFuzzyMatch = { String barcode, String seq, String qual, int from, int bcIndex, int slave -> // fuzzy match from seed match
    int nBadMMs = 0, nGoodMMs = 0
    int goodThreshold = maxGoodMMs[slave][bcIndex], badThreshold = maxBadMMs[slave][bcIndex]
    for (int i = 0; i < barcode.size(); i++) {
        if (!(redundancy2[(Character) barcode.charAt(i)].contains((Character) seq.charAt(from + i)))) {
            if (Util.qualFromSymbol(qual.charAt(i)) > lq)
                if (++nGoodMMs > goodThreshold)
                    break
                else if (++nBadMMs > badThreshold)
                    break
        }
    }
    nGoodMMs <= goodThreshold && nBadMMs <= badThreshold
}
def findMatch = { String barcode, Pattern seed, String seq, String qual, int bcIndex, int slave ->
    def seedOccurences = findAllMatches(seq, seed) // All seed occurences

    for (int i = 0; i < seedOccurences.size(); i++) {
        if (hasFuzzyMatch(barcode, seq, qual, seedOccurences[i], bcIndex, slave))  // Exhaustive search
            return seedOccurences[i] // till first best match
    }
    return -1
}

//========================
//  READ PROCESSING UTILS
//========================
// Overlap, select top quality nts for overlap zone
int maxConsMms = 2
double maxOverlapMismatchRatio = 0.1

def overlapReads = { String r1, String r2, String q1, String q2 ->
    String[] result = null

    for (int i = 0; i < maxOverlapOffset; i++) {
        if (i + overlapSeedSize > r2.length())
            return null

        def kmer = r2.substring(i, i + overlapSeedSize)
        def pattern = Pattern.compile(kmer)
        def matcher = pattern.matcher(r1)
        // Find last match
        int position
        while (matcher.find()) {
            position = matcher.start()
            if (position >= 0) {
                // Start fuzzy align
                boolean alignedAll = true
                int nConsMms = 0, nMms = 0, actualFuzzyOverlapSize = overlapFuzzySize

                for (int j = 0; j < overlapFuzzySize; j++) {
                    def posInR1 = position + overlapSeedSize + j, posInR2 = i + overlapSeedSize + j
                    if (posInR1 + 1 > r1.length() || posInR2 + 1 > r2.length()) {
                        actualFuzzyOverlapSize = j + 1
                        alignedAll = false
                        break     // went to end of r1
                    }
                    if (r1.charAt(posInR1) != r2.charAt(posInR2)) {
                        nMms++
                        if (++nConsMms >= maxConsMms)
                            break  // several consequent mismatches
                    } else {
                        nConsMms = 0 // zero counter
                    }
                }

                if (nConsMms < maxConsMms &&
                        (allowPartialOverlap || alignedAll) &&
                        (nMms / (double) actualFuzzyOverlapSize) <= maxOverlapMismatchRatio) {
                    // take best qual nts
                    def seq = new StringBuilder(r1.substring(0, position)),
                        qual = new StringBuilder(q1.substring(0, position))

                    int pos2 = i - 1
                    for (int j = position; j < r1.length(); j++) {
                        pos2++
                        if (pos2 == r2.length())
                            break // should not happen

                        seq.append(Util.qualFromSymbol(q1.charAt(j)) > Util.qualFromSymbol(q2.charAt(pos2)) ?
                                r1.charAt(j) : r2.charAt(pos2))

                        qual.append(Util.qualFromSymbol(q1.charAt(j)) > Util.qualFromSymbol(q2.charAt(pos2)) ?
                                q1.charAt(j) : q2.charAt(pos2))
                    }
                    for (int j = pos2 + 1; j < r2.length(); j++) {
                        // fill the remainder
                        seq.append(r2.charAt(j))

                        qual.append(q2.charAt(j))
                    }

                    // report overlap
                    result = new String[2]
                    result[0] = seq.toString()
                    result[1] = qual.toString()

                    return result
                }
            }
        }
    }

    result // failed
}

def flip = { String[] x, int i, int j ->
    String tmp
    tmp = x[i]
    x[i] = x[j]
    x[j] = tmp
}
def counters = new HashMap<String, AtomicLong[]>()
def readCounter = new AtomicLong(), goodReadCounter = new AtomicLong(), masterFirstCounter = new AtomicLong()
def overlapCounter = new AtomicLong()
def wrapRead = { String[] readData, StringBuilder[] umiData, int readIndex, String sampleId, String masterSampleId,
                 boolean good ->
    String umiHeader = umiData[0].size() > 0 ? " UMI:${umiData[0].toString()}:${umiData[1].toString()}" : ""
    readData[0] = readData[0] + umiHeader

    if (paired) {
        readData[3] = readData[3] + umiHeader

        if (readIndex > 0) { // flip if master is R2
            flip(readData, 0, 3)
            flip(readData, 1, 4)
            flip(readData, 2, 5)
        }

        if (rcReadMask[1]) { // rc slave if needed
            readData[4] = Util.revCompl(readData[4])
            readData[5] = readData[5].reverse()
        }
        readData[6] = sampleId
    } else readData[3] = sampleId

    if (rcReadMask[0]) {  // rc master if needed
        readData[1] = Util.revCompl(readData[1])
        readData[2] = readData[2].reverse()
    }

    // Record stats and try overlap
    long nReads = readCounter.incrementAndGet(),
         nGoodReads = good ? goodReadCounter.incrementAndGet() : goodReadCounter.get(),
         nMasterFirst = (good && (readIndex > 0)) ? masterFirstCounter.incrementAndGet() : masterFirstCounter.get(),
         overlapCount = overlapCounter.get()

    if (good) {
        counters.get(sampleId)[0].incrementAndGet()
        if (paired) {
            counters.get(sampleId)[1].incrementAndGet()

            if (overlap) {
                def overlapResult = overlapReads(readData[1], readData[4], readData[2], readData[5])

                if (overlapResult) {
                    overlapCounter.incrementAndGet()
                    readData[1] = overlapResult[0]
                    readData[2] = overlapResult[1]
                    readData[7] = "yes"

                    counters.get(sampleId)[2].incrementAndGet()
                }
            }
        }
    } else {
        counters.get(masterSampleId)[0].incrementAndGet() // where master was found
        if (paired)
            counters.get(sampleId)[1].incrementAndGet()
    }

    if (nReads % 25000 == 0)
        println "[${new Date()} $scriptName] Processed $nReads, " +
                "identified $nGoodReads (${((int) (10000 * (double) nGoodReads / (double) nReads)) / 100}%), " +
                (overlap ? "overlapped ${((int) (10000 * (double) overlapCount / (double) nGoodReads)) / 100}% of them, " : "") +
                "asymmetry (master first): $nMasterFirst (${((int) (10000 * (double) nMasterFirst / (double) nGoodReads)) / 100}%)"

    readData
}

//========================
//          BODY
//========================
def readQueue = new LinkedBlockingQueue<String[]>(4096)  // SINGLE: h1 r1 q1 ''       PAIRED: h1 r1 q1 h2 r2 q2 ''       ''
def writeQueue = new LinkedBlockingQueue<String[]>(4096) // SINGLE: h1 r1 q1 sampleId PAIRED: h1 r1 q1 h2 r2 q2 sampleId overlapped

if (!(new File(inputFileName1).exists())) {
    inputFileName1 = new File(barcodesFileName).parentFile.absolutePath + "/" + inputFileName1
}
if (inputFileName2 != BLANK_PATH && !(new File(inputFileName2).exists())) {
    inputFileName2 = new File(barcodesFileName).parentFile.absolutePath + "/" + inputFileName2
}

def reader1 = Util.getReader(inputFileName1), reader2 = (inputFileName2 == BLANK_PATH) ? null : Util.getReader(inputFileName2)
def writers = new HashMap<String, BufferedWriter[]>()
println "[${new Date()} $scriptName] Started processing for $inputFileName1, $inputFileName2"
new File(outputDir).mkdir()
[sampleIds, 'undef-s', 'undef-m'].flatten().each { String sampleId ->
    def writerTrio = new BufferedWriter[3]

    writerTrio[0] = Util.getWriter(paired ? "$outputDir/${sampleId}_R1.fastq" :
            "$outputDir/${sampleId}_R0.fastq", compressed, append)
    if (paired) {
        writerTrio[1] = Util.getWriter("$outputDir/${sampleId}_R2.fastq", compressed, append)
        if (overlap)
            writerTrio[2] = Util.getWriter("$outputDir/${sampleId}_R12.fastq", compressed, append)
    }

    writers.put(sampleId, writerTrio)

    def counterTrio = new AtomicLong[3]

    counterTrio[0] = new AtomicLong()
    if (paired) {
        counterTrio[1] = new AtomicLong()
        if (overlap)
            counterTrio[2] = new AtomicLong()
    }

    counters.put(sampleId, counterTrio)
}

int nProcessors = 2 * THREADS
def rDataSize = paired ? 8 : 4
int readsTaken = 0
def readThread = new Thread({  // Reading thread
    String header1
    while ((header1 = reader1.readLine()) != null) {
        def readData = new String[rDataSize]
        readData[0] = header1 + " R1"
        readData[1] = reader1.readLine()
        reader1.readLine()
        readData[2] = reader1.readLine()

        if (paired) {
            readData[3] = reader2.readLine() + " R2"
            readData[4] = reader2.readLine()
            reader2.readLine()
            readData[5] = reader2.readLine()
        }
        readQueue.put(readData)

        if (firstReadsToTake >= 0 && ++readsTaken > firstReadsToTake)
            break
    }
    for (int k = 0; k < nProcessors; k++)
        readQueue.put(new String[0]) // empty read - finish. tell ALL processors to stop
} as Runnable)
readThread.start()

def pool = Executors.newFixedThreadPool(THREADS)
def futures = new ArrayList<Future>()

def tgPattern = /^T{0,3}G{3,7}/, caPattern = /C{3,7}A{0,3}$/

for (int k = 0; k < nProcessors; k++) {
    futures.add(pool.submit(new Runnable() { // Processors
        @Override
        void run() {
            while (true) {
                def readData = readQueue.take()
                if (readData.length == 0)
                    break
                def umiData = new StringBuilder[2]
                umiData[0] = new StringBuilder()
                umiData[1] = new StringBuilder()

                // 1. Search master barcode
                int sampleIndex = -1, readIndex = -1, from = -1
                for (int i = 0; i < barcodeList[0].size(); i++) {
                    // R1
                    if ((from = findMatch(barcodeList[0][i], seeds[0][i],
                            readData[1], readData[2], i, 0)) >= 0) {
                        sampleIndex = i
                        readIndex = 0
                        break
                    }
                    // R2
                    if (!oriented && paired &&
                            ((from = findMatch(barcodeList[0][i], seeds[0][i],
                                    readData[4], readData[5], i, 0)) >= 0)) {
                        sampleIndex = i
                        readIndex = 3
                        break
                    }
                }

                if (from >= 0) {
                    // 1.2 Extract UMI from master
                    if (extractUMI) {
                        def seq = readData[readIndex + 1], qual = readData[readIndex + 2]
                        umiPositions[0][sampleIndex].each { int pos ->
                            umiData[0].append(seq.charAt(from + pos))
                            umiData[1].append(qual.charAt(from + pos))
                        }
                    }

                    // trim adapter
                    if (trimBc) {
                        int to = barcodeList[0][sampleIndex].size() + from
                        def seq = readData[readIndex + 1], qual = readData[readIndex + 2]
                        readData[readIndex + 1] = (from > trimSizeThreshold ? seq.substring(0, from) : "") +
                                ((seq.length() - to > trimSizeThreshold) ? seq.substring(to) : "")
                        readData[readIndex + 2] = (from > trimSizeThreshold ? qual.substring(0, from) : "") +
                                ((seq.length() - to > trimSizeThreshold) ? qual.substring(to) : "")
                        if (removeTS) {
                            seq = readData[readIndex + 1]
                            qual = readData[readIndex + 2]
                            def tsHit = (seq =~ tgPattern)
                            if (tsHit) {
                                readData[readIndex + 1] = seq.substring(tsHit.end())
                                readData[readIndex + 2] = qual.substring(tsHit.end())
                            }
                        }
                    }

                    // 2. Look for slave
                    if (paired) {
                        boolean slaveFound = barcodeList[1][sampleIndex] == 'N' // slave could be skipped
                        if (!slaveFound) {
                            // Look in other read
                            from = findMatch(barcodeList[1][sampleIndex], seeds[1][sampleIndex],
                                    readData[4 - readIndex], readData[5 - readIndex], sampleIndex, 1)

                            slaveFound = from >= 0

                            if (slaveFound) {
                                def seq = readData[4 - readIndex], qual = readData[5 - readIndex]

                                // 2.1 Extract UMI from slave
                                if (extractUMI) {
                                    umiPositions[1][sampleIndex].each { int pos ->
                                        umiData[0].append(seq.charAt(from + pos))
                                        umiData[1].append(qual.charAt(from + pos))
                                    }
                                }

                                // trim adapter
                                if (trimBc) {
                                    int to = barcodeList[1][sampleIndex].size() + from
                                    readData[4 - readIndex] = (from > trimSizeThreshold ? seq.substring(0, from) : "") +
                                            ((seq.length() - to > trimSizeThreshold) ? seq.substring(to) : "")
                                    readData[5 - readIndex] = (from > trimSizeThreshold ? qual.substring(0, from) : "") +
                                            ((seq.length() - to > trimSizeThreshold) ? qual.substring(to) : "")

                                    if (removeTS) {
                                        seq = readData[4 - readIndex]
                                        qual = readData[5 - readIndex]
                                        def tsHit = (seq =~ caPattern)
                                        if (tsHit) {
                                            readData[readIndex + 1] = seq.substring(0, tsHit.start())
                                            readData[readIndex + 2] = qual.substring(0, tsHit.start())
                                        }
                                    }
                                }
                            }
                        }
                        if (slaveFound) {
                            // 2.2 Report full match
                            writeQueue.put(wrapRead(readData, umiData, readIndex,
                                    sampleIds[sampleIndex], sampleIds[sampleIndex], true))
                        } else {
                            // 2.3 Report failed slave
                            def wrappedRead = wrapRead(readData, umiData, readIndex,
                                    'undef-s', sampleIds[sampleIndex], false)
                            if (!noUndef)
                                writeQueue.put(wrappedRead)
                        }
                    } else { // Single-end, master match
                        writeQueue.put(wrapRead(readData, umiData, readIndex,
                                sampleIds[sampleIndex], sampleIds[sampleIndex], true))
                    }
                } else {
                    // 1.3 Report failed master
                    def wrappedRead = wrapRead(readData, umiData, readIndex,
                            'undef-m', 'undef-m', false)
                    if (!noUndef)
                        writeQueue.put(wrappedRead)
                }
            }
        }
    } as Runnable))
}

def writeThread = new Thread({  // Writing thread
    String[] result
    while (true) {
        result = writeQueue.take()

        if (result.length == 0)
            break

        def sampleId = paired ? result[6] : result[3]

        def writerTrio = writers.get(sampleId)

        if (!paired)
            writerTrio[0].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
        else {
            if (overlap && result[7]) {
                writerTrio[2].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
            } else {
                writerTrio[0].writeLine(result[0] + "\n" + result[1] + "\n+\n" + result[2])
                writerTrio[1].writeLine(result[3] + "\n" + result[4] + "\n+\n" + result[5])
            }
        }
    }
} as Runnable)
writeThread.start()

readThread.join() // wait for read to finish

futures.each { it.get() } // wait for processing to finish
writeQueue.put(new String[0]) // tell writers this is last one
pool.shutdown()

writeThread.join()  // wait for write to finish

writers.values().each { // don't forget to flush
    it[0].close()
    if (paired) {
        it[1].close()
        if (overlap)
            it[2].close()
    }
}

def fastqPrefix = compressed ? "fastq.gz" : "fastq"

def filelistFile = new File("$outputDir/checkout.filelist.txt")
def filelistExists = filelistFile.exists()
def filelistWriter = filelistFile.newWriter(append)

def existingSampleIdList = new HashSet<String>()
if (!filelistExists) {
    filelistWriter.println("#SAMPLE_ID\tSAMPLE_TYPE\tCHECKOUT_FASTQ1\tCHECKOUT_FASTQ2")
} else if (append) {
    new File("$outputDir/checkout.filelist.txt").splitEachLine("\t") {
        existingSampleIdList.add(it[0])
    }
}

new HashSet(sampleIds).each { String sampleId -> // only unique
    if (!existingSampleIdList.contains(sampleId)) {
        filelistWriter.println(sampleId +
                (paired ? "\tpaired\t" : "\tunpaired\t") +
                (paired ? new File("$outputDir/${sampleId}_R1.$fastqPrefix").absolutePath :
                        new File("$outputDir/${sampleId}_R0.$fastqPrefix").absolutePath)
                + "\t" +
                (paired ? new File("$outputDir/${sampleId}_R2.$fastqPrefix").absolutePath : BLANK_PATH))
        if (overlap)
            filelistWriter.println(sampleId + "\toverlapped\t" +
                    new File("$outputDir/${sampleId}_R12.$fastqPrefix").absolutePath + "\t$BLANK_PATH")
    }
}

filelistWriter.close()

def logFile = new File("$outputDir/checkout.log.txt")
def logExists = logFile.exists()
def logWriter = logFile.newWriter(append)
if (!logExists) {
    logWriter.println("#INPUT_FILE_1\tINPUT_FILE_2\tSAMPLE\tMASTER\tMASTER+SLAVE\tOVERLAPPED")
}

counters.each {
    logWriter.println(inputFileName1 + "\t" + inputFileName2 + "\t" + it.key + "\t" + it.value.collect().join("\t"))
}

logWriter.close()

Util.printCmd(outputDir + "/checkout.cmd.txt")

println "[${new Date()} $scriptName] Finished"
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
