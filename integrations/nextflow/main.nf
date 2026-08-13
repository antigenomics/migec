#!/usr/bin/env nextflow
// A runnable entry point, so the modules can be exercised without writing a pipeline first.
//
//   nextflow run integrations/nextflow -profile docker \
//       --mode ctdna --input 'data/*_R{1,2}.fq.gz' --preset tso500 \
//       --fasta ref.fa --outdir results/
//
//   nextflow run integrations/nextflow --mode airr --input 'data/*_R{1,2}.fq.gz' --preset migec
//   nextflow run integrations/nextflow --mode consensus --input '...' --preset 10x-v2
//
// `-stub-run` walks the whole graph with nothing installed.

nextflow.enable.dsl = 2

include { MIGEC       } from './subworkflows/migec/main'
include { MIGEC_CTDNA } from './subworkflows/migec_ctdna/main'
include { MIGEC_AIRR  } from './subworkflows/migec_airr/main'

params.mode      = 'consensus'      // consensus | ctdna | airr
params.input     = null
params.outdir    = 'results'
params.fasta     = null
params.preset    = null
params.single_end = false

workflow {
    if (!params.input) {
        error "--input is required, e.g. --input 'data/*_R{1,2}.fq.gz'"
    }

    // A preset given on the command line becomes the per-sample default. Per-sample keys in a
    // samplesheet still win, which is what lets one run mix chemistries.
    ch_reads = params.single_end
        ? Channel.fromPath(params.input).map { f -> [ [ id: f.simpleName, preset: params.preset,
                                                        single_end: true ], f ] }
        : Channel.fromFilePairs(params.input).map { id, files -> [ [ id: id,
                                                                     preset: params.preset ], files ] }

    ch_whitelist = params.getOrDefault('migec_cell_whitelist', null)
        ? Channel.value(file(params.migec_cell_whitelist))
        : Channel.value(file("$projectDir/assets/NO_FILE", checkIfExists: false))

    if (params.mode == 'ctdna') {
        if (!params.fasta) {
            error "--mode ctdna needs --fasta (and a .fai beside it)"
        }
        ch_fasta = Channel.value(file(params.fasta))
        ch_fai   = Channel.value(file("${params.fasta}.fai"))
        MIGEC_CTDNA(ch_reads, ch_fasta, ch_fai, ch_whitelist)
        MIGEC_CTDNA.out.vcf.view { meta, vcf -> "${meta.id}: ${vcf}" }
    }
    else if (params.mode == 'airr') {
        MIGEC_AIRR(ch_reads, ch_whitelist)
        MIGEC_AIRR.out.airr.view { meta, airr -> "${meta.id}: ${airr}" }
    }
    else if (params.mode == 'consensus') {
        MIGEC(ch_reads, ch_whitelist)
        MIGEC.out.consensus.view { meta, fq -> "${meta.id}: ${fq}" }
    }
    else {
        error "--mode must be consensus, ctdna or airr; got '${params.mode}'"
    }
}
