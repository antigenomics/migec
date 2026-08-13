// Immune repertoire, end to end: reads -> molecules -> AIRR clonotypes.
//
// The deliverable here is a COUNT, not a genotype, which changes what the barcode is for. Every
// PCR artefact is a valid-looking rare clonotype and there is no reference to contradict it, so
// counting reads counts the artefacts; counting molecules does not. `duplicate_count` in the AIRR
// output is a molecule count, and that is what a diversity estimate should use.
//
// Note: this drops straight into nf-core/airrflow in place of its own presto/UMI steps, because
// what it emits is ordinary AIRR TSV.

include { MIGEC           } from '../migec/main'
include { DOWNSTREAM_ARDA } from '../../modules/downstream/arda/main'

workflow MIGEC_AIRR {
    take:
    ch_reads      // channel: [ val(meta), [ reads ] ]
    ch_whitelist  // value channel: a cell whitelist, or file('NO_FILE')

    main:
    ch_versions = Channel.empty()

    MIGEC(ch_reads, ch_whitelist)
    ch_versions = ch_versions.mix(MIGEC.out.versions)

    DOWNSTREAM_ARDA(MIGEC.out.consensus)
    ch_versions = ch_versions.mix(DOWNSTREAM_ARDA.out.versions.first())

    emit:
    airr      = DOWNSTREAM_ARDA.out.airr
    consensus = MIGEC.out.consensus
    molecules = MIGEC.out.molecules
    cells     = MIGEC.out.cells
    qc        = MIGEC.out.qc
    versions  = ch_versions
}
