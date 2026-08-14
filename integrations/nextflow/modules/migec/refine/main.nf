// migec refine: correct barcode errors, call cells, hand the reads on.
//
// MEMORY scales with DISTINCT BARCODES, never with reads: refine holds the barcode TABLE and
// streams the reads past it. ~96 B per distinct barcode, so ~10 GB at 10^8 distinct barcodes and
// nothing at all proportional to file size. Budget by expected molecules, not by FASTQ bytes --
// a 200 GB shallow run and a 200 GB deep one need wildly different memory here.
//
// CPU: the neighbourhood scan and the read rewrite both thread, and the output is byte-identical
// at any `--threads`, so `task.cpus` is free to vary between attempts.

process MIGEC_REFINE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_high_memory'   // barcodes, not reads -- see above

    conda "${moduleDir}/../environment.yml"
    // One place, so a release bumps one line rather than four. Override with
    // `--migec_container` when you build your own image.
    container params.getOrDefault('migec_container', 'migec:2.3.0')

    input:
    tuple val(meta), path(reads)
    path cell_whitelist   // optional; stage an empty file or [] when there is none

    output:
    tuple val(meta), path("ref/${meta.id}.fq.gz"),      emit: reads
    tuple val(meta), path("ref/${meta.id}.barcodes.tsv"), emit: barcodes
    tuple val(meta), path("ref/${meta.id}.bins.tsv"),   emit: bins
    tuple val(meta), path("ref/${meta.id}.rank.tsv"),   emit: rank
    tuple val(meta), path("ref/${meta.id}.cells.tsv"),  emit: cells, optional: true
    tuple val(meta), path("ref/refine.*.tsv"),          emit: qc
    tuple val(meta), path("ref/refine.json"),           emit: report
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cells  = meta.expect_cells ?: params.getOrDefault('migec_expect_cells', 3000)
    def wl     = cell_whitelist && cell_whitelist.name != 'NO_FILE'
                     ? "--cell-whitelist ${cell_whitelist}" : ''
    def limit  = params.getOrDefault('migec_limit_read', 0)
    def limitu = params.getOrDefault('migec_limit_umi', 0)
    """
    migec refine ${reads} \\
        --sample ${prefix} \\
        --expect-cells ${cells} \\
        ${wl} \\
        ${limit ? "--limit-read ${limit}" : ''} \\
        ${limitu ? "--limit-umi ${limitu}" : ''} \\
        --threads ${task.cpus} \\
        ${params.getOrDefault('migec_refine_args', '')} \\
        -o ref/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info | awk '/^migec  /{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ref
    echo | gzip > ref/${prefix}.fq.gz
    for f in ${prefix}.barcodes.tsv ${prefix}.bins.tsv ${prefix}.rank.tsv refine.coverage.tsv; do
        touch ref/\$f
    done
    echo '{}' > ref/refine.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info 2>/dev/null | awk '/^migec  /{print \$2}' || echo unknown)
    END_VERSIONS
    """
}
