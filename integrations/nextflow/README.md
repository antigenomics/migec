# migec as an nf-core-style local module

```
modules/migec/checkout/main.nf     reads     -> tagged FASTQ + QC tables
modules/migec/refine/main.nf       tagged    -> corrected FASTQ + barcode table + cell calls
modules/migec/assemble/main.nf     corrected -> one consensus per molecule
modules/migec/plot/main.nf         the tables -> SVG figures (needs gnuplot)
subworkflows/migec/main.nf         the three chained, with the payload mate handled
nextflow.config                    defaults, every one read with getOrDefault
```

```groovy
include { MIGEC } from './integrations/nextflow/subworkflows/migec/main'

workflow {
    ch_reads = Channel.fromFilePairs(params.input)
        .map { id, files -> [ [ id: id, preset: '10x-v2', payload_mate: 2 ], files ] }

    MIGEC(ch_reads, file(params.cell_whitelist ?: 'NO_FILE'))

    MIGEC.out.consensus.view()
}
```

Per-sample keys in `meta` win over the `params.migec_*` defaults, so one run can mix chemistries:
`bc_pattern`, `preset`, `read_structure`, `read_structure2`, `barcodes`, `max_offset`,
`payload_mate`, `expect_cells`, `rt_error`, `contig`, `fast`.

**The full write-up is [`docs/nextflow.rst`](../../docs/nextflow.rst)** — why it is three processes
and not one, the resource labels, and the three things that go wrong (the payload mate, a
per-sample `false`, and `--rt-error` naming a chemistry rather than a number). It is one file so
the throughput numbers and the caveats cannot drift apart from the rest of the docs.

`nextflow run -stub-run` walks the whole graph without migec installed; every process has a `stub:`
block. Note that neither the stub run nor a real run has been executed here — nextflow is not
installed on the machine this was written on, so treat your first pipeline run as the verification.
