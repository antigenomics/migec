// migec assemble: collapse each molecule's reads into one consensus.
//
// MEMORY is bounded by the largest resident BUCKET, times the number of workers -- the input is
// range-partitioned into at least 16 buckets and each worker sorts one at a time. It does not
// scale with the library. CPU threads over the buckets, and the output is byte-identical at any
// `--threads` because the bucket count is a property of the input, never of `task.cpus`.
//
// `--rt-error` names the chemistry rather than taking a number: `rt` (1e-4, caps at Q40, anything
// with a reverse transcription step), `medium` (1e-5), `high` (1e-6, a proofreading polymerase and
// no RT). It is the ONE-MOLECULE floor -- 10x's Q60 is for bases covered by two or more UMIs, and
// combining molecules is arda's job, downstream of this.

process MIGEC_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container "migec:2.0.0a3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("asm/${meta.id}.consensus.fq.gz"), emit: consensus
    tuple val(meta), path("asm/${meta.id}.mig.tsv"),         emit: molecules
    tuple val(meta), path("asm/assemble.coverage.tsv"),      emit: qc
    tuple val(meta), path("asm/assemble.json"),              emit: report
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rt     = meta.rt_error ?: params.getOrDefault('migec_rt_error', 'rt')
    def contig = (meta.contig ?: params.getOrDefault('migec_contig', false)) ? '--contig' : ''
    def fast   = (meta.fast ?: params.getOrDefault('migec_fast', false)) ? '--fast' : ''
    """
    migec assemble ${reads} \\
        --sample ${prefix} \\
        --rt-error ${rt} \\
        ${contig} ${fast} \\
        --threads ${task.cpus} \\
        ${params.getOrDefault('migec_assemble_args', '')} \\
        -o asm/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info | awk '/^migec  /{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p asm
    echo | gzip > asm/${prefix}.consensus.fq.gz
    touch asm/${prefix}.mig.tsv asm/assemble.coverage.tsv
    echo '{}' > asm/assemble.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info 2>/dev/null | awk '/^migec  /{print \$2}' || echo unknown)
    END_VERSIONS
    """
}
