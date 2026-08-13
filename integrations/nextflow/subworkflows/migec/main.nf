// checkout -> refine -> assemble, wired. This is what a pipeline includes.
//
// Never: THE PAYLOAD MATE IS A DECLARATION, NOT A GUESS. checkout writes <sample>_R1/_R2 when the
// input was paired, and on a droplet chemistry the barcode read carries no cDNA at all -- so which
// mate the later stages see has to be said. `params.migec_payload_mate` (or meta.payload_mate)
// selects it and defaults to 1, because a bulk amplicon has only one.

include { MIGEC_CHECKOUT } from '../../modules/migec/checkout/main'
include { MIGEC_REFINE   } from '../../modules/migec/refine/main'
include { MIGEC_ASSEMBLE } from '../../modules/migec/assemble/main'
include { MIGEC_PLOT     } from '../../modules/migec/plot/main'

workflow MIGEC {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_whitelist      // value channel: a cell whitelist, or file('NO_FILE')

    main:
    ch_versions = Channel.empty()

    MIGEC_CHECKOUT(ch_reads)
    ch_versions = ch_versions.mix(MIGEC_CHECKOUT.out.versions.first())

    // One file when single-end, two when paired: pick the mate that carries sequence.
    ch_payload = MIGEC_CHECKOUT.out.reads.map { meta, files ->
        def list = files instanceof List ? files : [files]
        def mate = (meta.payload_mate ?: params.getOrDefault('migec_payload_mate', 1)) as int
        def pick = list.size() > 1 ? list.find { it.name.endsWith("_R${mate}.fq.gz") } : list[0]
        if (pick == null) {
            error "migec: no _R${mate} output for ${meta.id} -- checkout wrote ${list*.name}"
        }
        [ meta, pick ]
    }

    MIGEC_REFINE(ch_payload, ch_whitelist)
    ch_versions = ch_versions.mix(MIGEC_REFINE.out.versions.first())

    MIGEC_ASSEMBLE(MIGEC_REFINE.out.reads)
    ch_versions = ch_versions.mix(MIGEC_ASSEMBLE.out.versions.first())

    // Every stage's tables, per sample, drawn in one pass.
    ch_tables = MIGEC_CHECKOUT.out.qc
        .mix(MIGEC_REFINE.out.qc, MIGEC_REFINE.out.rank,
             MIGEC_ASSEMBLE.out.qc, MIGEC_ASSEMBLE.out.molecules)
        .groupTuple()
        .map { meta, files -> [ meta, files.flatten() ] }

    if (params.getOrDefault('migec_plots', true)) {
        MIGEC_PLOT(ch_tables)
        ch_versions = ch_versions.mix(MIGEC_PLOT.out.versions.first())
    }

    emit:
    consensus = MIGEC_ASSEMBLE.out.consensus   // [ meta, *.consensus.fq.gz ]
    molecules = MIGEC_ASSEMBLE.out.molecules   // [ meta, *.mig.tsv ]
    barcodes  = MIGEC_REFINE.out.barcodes
    cells     = MIGEC_REFINE.out.cells
    qc        = ch_tables
    versions  = ch_versions
}
