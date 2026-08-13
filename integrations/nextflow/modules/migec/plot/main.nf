// migec plot: QC figures from the tables the stages already wrote.
//
// Its own process because it needs gnuplot and nothing else does, so a pipeline without gnuplot
// drops this one module and loses the figures rather than the run. It reads no reads: give it a
// directory of TSVs and it writes SVGs beside them, so it is cheap, resumable, and never on the
// critical path.

process MIGEC_PLOT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/../environment.yml"
    // One place, so a release bumps one line rather than four. Override with
    // `--migec_container` when you build your own image.
    container params.getOrDefault('migec_container', 'migec:2.1.0')

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("plots/*.svg"), emit: figures, optional: true
    tuple val(meta), path("plots/*.gp"),  emit: scripts
    path "versions.yml",                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p tables && cp ${tables} tables/ 2>/dev/null || true
    migec plot tables/ -o plots/ --format ${params.getOrDefault('migec_plot_format', 'svg')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info | awk '/^migec  /{print \$2}')
        gnuplot: \$(gnuplot --version 2>/dev/null | awk '{print \$2}' || echo "not installed")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p plots && touch plots/coverage.gp
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        migec: \$(migec info 2>/dev/null | awk '/^migec  /{print \$2}' || echo unknown)
    END_VERSIONS
    """
}
