// Align the consensus, carrying the molecule tags into the BAM.
//
// Never: THE COMMENT FLAG IS ALIGNER-SPECIFIC AND THE WRONG ONE IS SILENT-ADJACENT.
// minimap2 and `minibwa map` spell it `-y`; bwa's `mem` spells it `-C`. Each rejects the other
// with a non-zero exit rather than dropping the tags quietly, which is the only reason this is
// safe to parameterise -- see docs/downstream.rst, where all three were run.
//
// One consensus record is one molecule, so this aligns ONE record per molecule rather than one per
// read. On a library at 15 reads/molecule that is a fifteenth of the alignment work, and the
// record being aligned is the error-corrected one.

process DOWNSTREAM_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.getOrDefault('align_container',
                                  'quay.io/biocontainers/minimap2:2.28--he4a0461_0')

    input:
    tuple val(meta), path(consensus)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "versions.yml",                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def aligner = meta.aligner ?: params.getOrDefault('downstream_aligner', 'minimap2')
    def cmd
    if (aligner == 'minimap2') {
        cmd = "minimap2 -ax sr -y -t ${task.cpus} ${fasta} ${consensus}"
    }
    else if (aligner == 'bwa') {
        cmd = "bwa mem -C -t ${task.cpus} ${fasta} ${consensus}"
    }
    else if (aligner == 'minibwa') {
        cmd = "minibwa map -y -t ${task.cpus} ${fasta} ${consensus}"
    }
    else {
        error "downstream_aligner must be minimap2, bwa or minibwa; got '${aligner}'"
    }
    """
    ${cmd} | samtools sort -@ ${task.cpus} -o ${prefix}.bam -
    samtools index ${prefix}.bam

    # The tags are the point of this step, so fail here rather than three processes later.
    if ! samtools view ${prefix}.bam | head -1 | grep -q 'MI:Z:'; then
        echo "ERROR: no MI:Z: tag in ${prefix}.bam -- the aligner dropped the FASTQ comment." >&2
        echo "minimap2/minibwa need -y, bwa mem needs -C." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aligner: ${aligner}
        samtools: \$(samtools --version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aligner: stub
    END_VERSIONS
    """
}
