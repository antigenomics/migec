// ctDNA / rare somatic variants, end to end: reads -> molecules -> BAM -> VCF.
//
// The shape this encodes, from docs/variants.rst:
//
//   collapse FIRST, then align once, then call with a STANDARD caller.
//
// Aligning raw reads and grouping on (position, UMI) is the other order in use -- fgbio,
// UMI-tools, UMIErrorCorrect -- and it costs one alignment per READ instead of one per MOLECULE,
// with the aligner seeing uncorrected sequence. The position is worth real key bits when the
// barcode is short; assemble recovers that discriminating power from the payload instead, by
// sub-clustering a barcode group whose reads carry two co-segregating haplotypes.
//
// Never: DO NOT ADD A DEDUPLICATION STEP HERE. There is exactly one record per molecule after
// MIGEC. MarkDuplicates, UMI-VarCal or UMIErrorCorrect placed after this delete real molecules.

include { MIGEC                   } from '../migec/main'
include { DOWNSTREAM_ALIGN        } from '../../modules/downstream/align/main'
include { DOWNSTREAM_CALLVARIANTS } from '../../modules/downstream/callvariants/main'

workflow MIGEC_CTDNA {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_fasta      // value channel: reference FASTA
    ch_fai        // value channel: its .fai
    ch_whitelist  // value channel: a cell whitelist, or file('NO_FILE')

    main:
    ch_versions = Channel.empty()

    MIGEC(ch_reads, ch_whitelist)
    ch_versions = ch_versions.mix(MIGEC.out.versions)

    DOWNSTREAM_ALIGN(MIGEC.out.consensus, ch_fasta)
    ch_versions = ch_versions.mix(DOWNSTREAM_ALIGN.out.versions.first())

    DOWNSTREAM_CALLVARIANTS(DOWNSTREAM_ALIGN.out.bam, ch_fasta, ch_fai)
    ch_versions = ch_versions.mix(DOWNSTREAM_CALLVARIANTS.out.versions.first())

    emit:
    vcf       = DOWNSTREAM_CALLVARIANTS.out.vcf
    bam       = DOWNSTREAM_ALIGN.out.bam
    consensus = MIGEC.out.consensus
    molecules = MIGEC.out.molecules
    qc        = MIGEC.out.qc
    versions  = ch_versions
}
