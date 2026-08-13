// Assign V/D/J and call clonotypes from the consensus, emitting AIRR.
//
// Note: THE AIRR `sequence_id` IS THE MOLECULE ID. arda reads the consensus through `dnaio`,
// which drops FASTQ comments -- which is exactly why `assemble` writes the identity into the read
// NAME (`<sample>.<cell>.<umi>`) as well as into SAM tags. The name stands alone, so nothing is
// lost here even though the tags are.
//
// Note: `duplicate_count` in the output is therefore a MOLECULE count, not a read count. That is
// the number a diversity estimate should be computed from; a read count is an amplification
// artefact wearing an abundance's clothes.

process DOWNSTREAM_ARDA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.getOrDefault('arda_container', 'arda:latest')

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
        --species ${species} \\
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
