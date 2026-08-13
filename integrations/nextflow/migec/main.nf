// migec: UMI-tagged reads -> one consensus per molecule, as a drop-in nf-core-style local module.
// One `checkout` -> `refine` -> `assemble` chain per sample; publishes to ${params.outdir}/migec/.
// See ./README.md for wiring into nf-core/airrflow (or any pipeline that hands you FASTQ pairs).
//
// Never: THE BARCODE READ IS NOT ALWAYS THE PAYLOAD READ. On 10x the barcode read is 26 nt of cell
// barcode and UMI and nothing else, and the whole cDNA is on the other mate -- so `checkout` is
// given both, and `refine`/`assemble` then run on the mate that carries sequence. This module
// picks that mate from `params.migec_payload_mate` rather than assuming R1, because assuming R1
// silently produces empty consensuses on every droplet chemistry.

process MIGEC {
    tag "$meta.id"

    // migec is CPU-bound in `checkout` (1.2 M reads/s at 16 threads on an M-series laptop) and
    // single-threaded afterwards. `refine` streams the reads three times and holds the barcode
    // TABLE -- ~96 B per distinct barcode, so ~10 GB at 10^8 distinct barcodes and nothing at all
    // proportional to read count. `assemble` holds one range-partitioned bucket at a time.
    //
    // MEMORY scales with DISTINCT BARCODES, never with reads. Measured: 3.16 M 10x VDJ reads over
    // 311 k barcodes peaked 658 MB; the same read count over 4x the barcodes peaks ~4x the table.
    // Budget by expected molecules, not by FASTQ size.
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "migec:2.0.0a1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.consensus.fq.gz"),   emit: consensus
    tuple val(meta), path("*.mig.tsv"),           emit: molecules
    tuple val(meta), path("*.barcodes.tsv"),      emit: barcodes
    tuple val(meta), path("*.cells.tsv"),         emit: cells, optional: true
    tuple val(meta), path("*.rank.tsv"),          emit: rank
    tuple val(meta), path("*.bins.tsv"),          emit: bins
    tuple val(meta), path("checkout.*.tsv"),      emit: checkout_qc
    tuple val(meta), path("*.json"),              emit: report
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "${meta.id}"
    // getOrDefault throughout, never a bare params read. A module must run when its own
    // nextflow.config was not includeConfig'd -- a bare read emits "Access to undefined
    // parameter" on every one, which is a WARN normally and a hard failure under strict mode.
    def pattern  = meta.bc_pattern ?: params.getOrDefault('migec_bc_pattern', null)
    def rs       = meta.read_structure ?: params.getOrDefault('migec_read_structure', null)
    def rs2      = meta.read_structure2 ?: params.getOrDefault('migec_read_structure2', null)
    def layout   = rs ? "--read-structure '${rs}'" + (rs2 ? " --read-structure2 '${rs2}'" : '')
                      : "--bc-pattern '${pattern}'"
    def offset   = meta.max_offset != null ? meta.max_offset
                                           : params.getOrDefault('migec_max_offset', -1)
    def payload  = meta.payload_mate ?: params.getOrDefault('migec_payload_mate', 1)
    def cells    = params.getOrDefault('migec_expect_cells', 3000)
    def rt       = params.getOrDefault('migec_rt_error', 1e-4)
    def r1       = reads instanceof List ? reads[0] : reads
    def r2       = (reads instanceof List && reads.size() > 1) ? reads[1] : ''
    def whitelist = params.getOrDefault('migec_cell_whitelist', null)
    def wl       = whitelist ? "--cell-whitelist ${whitelist}" : ''
    // `checkout` writes <sample>_R1/_R2 when paired and <sample> when not.
    def stage    = r2 ? "co/${prefix}_R${payload}.fq.gz" : "co/${prefix}.fq.gz"
    """
    migec checkout ${r1} ${r2} \\
        ${layout} \\
        --sample ${prefix} \\
        --max-offset ${offset} \\
        --threads ${task.cpus} \\
        ${params.getOrDefault('migec_checkout_args', '')} \\
        -o co/

    migec refine ${stage} \\
        --sample ${prefix} \\
        --expect-cells ${cells} \\
        ${wl} \\
        ${params.getOrDefault('migec_refine_args', '')} \\
        -o ref/

    migec assemble ref/${prefix}.fq.gz \\
        --sample ${prefix} \\
        --rt-error ${rt} \\
        ${params.getOrDefault('migec_assemble_args', '')} \\
        -o asm/

    # Flat, so the emit globs above do not have to know the stage layout.
    mv co/checkout.*.tsv co/*.json . 2>/dev/null || true
    mv ref/*.tsv ref/*.json asm/*.tsv asm/*.json asm/*.consensus.fq.gz . 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info | awk '/^migec  /{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.consensus.fq.gz
    for f in ${prefix}.mig.tsv ${prefix}.barcodes.tsv ${prefix}.rank.tsv ${prefix}.bins.tsv \\
             checkout.summary.tsv; do touch \$f; done
    touch ${prefix}.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info 2>/dev/null | awk '/^migec  /{print \$2}' || echo unknown)
    END_VERSIONS
    """
}
