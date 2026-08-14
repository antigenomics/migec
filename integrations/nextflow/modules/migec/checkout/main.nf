// migec checkout: find the barcode, extract it, trim, demultiplex. One process, nf-core style.
//
// Split out from the old single MIGEC process because the three stages have different shapes and
// a pipeline needs to say so: checkout is CPU-bound over READS, refine is bound by DISTINCT
// BARCODES, assemble by the largest bucket. One process meant one resource label for all three,
// one retry for all three, and no resume between them -- a failed assemble re-ran the whole
// demultiplex.
//
// Never: THE BARCODE READ IS NOT ALWAYS THE PAYLOAD READ. On 10x, R1 is 26 nt of cell barcode and
// UMI and nothing else, and the cDNA is entirely on R2 -- so checkout is given both mates and the
// later stages run on the one that carries sequence. `payload_mate` says which; assuming R1
// silently produces empty consensuses on every droplet chemistry.

process MIGEC_CHECKOUT {
    tag "$meta.id"
    label 'process_high'   // the only stage whose cost scales with reads rather than barcodes

    conda "${moduleDir}/../environment.yml"
    // One place, so a release bumps one line rather than four. Override with
    // `--migec_container` when you build your own image.
    container params.getOrDefault('migec_container', 'migec:2.3.0')

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("co/${meta.id}*.fq.gz"), emit: reads
    tuple val(meta), path("co/checkout.*.tsv"),    emit: qc
    tuple val(meta), path("co/checkout.json"),     emit: report
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // getOrDefault throughout, never a bare params read: a module has to run when its own
    // nextflow.config was not includeConfig'd, and a bare read of an undefined parameter is a
    // WARN normally and a hard failure under strict mode.
    def pattern = meta.bc_pattern ?: params.getOrDefault('migec_bc_pattern', null)
    def preset  = meta.preset ?: params.getOrDefault('migec_preset', null)
    def rs      = meta.read_structure ?: params.getOrDefault('migec_read_structure', null)
    def rs2     = meta.read_structure2 ?: params.getOrDefault('migec_read_structure2', null)
    def sheet   = meta.barcodes ?: params.getOrDefault('migec_barcodes', null)
    def layout  = sheet   ? "-b ${sheet}"
                : preset  ? "--preset '${preset}'"
                : rs      ? "--read-structure '${rs}'" + (rs2 ? " --read-structure2 '${rs2}'" : '')
                          : "--bc-pattern '${pattern}'"
    // Passed only when set. The default is automatic -- a caret, a slice list, a read structure or
    // a pattern with nothing to score all anchor at the first base -- and passing -1 anyway
    // reinstates the refusal the anchor exists to avoid.
    def maxoff  = meta.max_offset != null ? meta.max_offset
                                          : params.getOrDefault('migec_max_offset', null)
    def offset  = maxoff != null ? "--max-offset ${maxoff}" : ''
    def limit   = params.getOrDefault('migec_limit_read', 0)
    def limited = limit ? "--limit-read ${limit}" : ''
    def r1      = reads instanceof List ? reads[0] : reads
    def r2      = (reads instanceof List && reads.size() > 1) ? reads[1] : ''
    """
    migec checkout ${r1} ${r2} \\
        ${layout} \\
        --sample ${prefix} \\
        ${offset} \\
        ${limited} \\
        --threads ${task.cpus} \\
        ${params.getOrDefault('migec_checkout_args', '')} \\
        -o co/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info | awk '/^migec  /{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p co
    echo | gzip > co/${prefix}.fq.gz
    for f in summary coverage umi_composition barcode_space umi_quality trimming; do
        touch co/checkout.\$f.tsv
    done
    echo '{}' > co/checkout.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info 2>/dev/null | awk '/^migec  /{print \$2}' || echo unknown)
    END_VERSIONS
    """
}
