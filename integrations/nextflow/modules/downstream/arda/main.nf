// Assign V/D/J and call clonotypes from the consensus, emitting AIRR.
//
// AIRR `sequence_id` identifies an emitted consensus RECORD. The read name carries
// <sample>[.<cell>].<umi>[.c<contig>][.<molecule>] even when comments are discarded.
// Separate paired contigs share the UMI/molecule branch and must not be counted twice.
// `duplicate_count` is native record abundance, not guaranteed physical molecule abundance.
// Molecular summaries must group qualified UMI/collision identities across contig records.

process DOWNSTREAM_ARDA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.getOrDefault('arda_container', null)

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("${meta.id}*.tsv"), emit: airr
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def species = meta.species ?: params.getOrDefault('arda_species', 'human')
    def args    = params.getOrDefault('arda_args', '')
    """
    arda amplicon \\
        --r1 ${consensus} \\
        --organism ${species} \\
        --threads ${task.cpus} \\
        ${args} \\
        -p ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: \$(arda --version 2>&1 | tail -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.airr.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arda: stub
    END_VERSIONS
    """
}
