# migec as an nf-core-style local module

Three processes and a subworkflow, laid out the way nf-core expects:

```
modules/migec/checkout/main.nf     reads  -> tagged FASTQ + QC tables
modules/migec/refine/main.nf       tagged -> corrected FASTQ + barcode table + cell calls
modules/migec/assemble/main.nf     corrected -> one consensus per molecule
modules/migec/plot/main.nf         the tables -> SVG figures (needs gnuplot)
subworkflows/migec/main.nf         the three chained, with the payload mate handled
nextflow.config                    defaults, all read with getOrDefault
```

## Why three processes and not one

The stages have different shapes, and one process forces one answer for all three:

| stage | scales with | threads | memory |
|---|---|---|---|
| `checkout` | reads | yes, 1.06 M reads/s at 16 | chunk-bounded, plus the UMI counters |
| `refine` | **distinct barcodes** | yes, 979 k reads/s at 16 | the barcode table, ~96 B each |
| `assemble` | reads, then buckets | yes, 1.49 M reads/s at 16 | one bucket per worker |

Splitting them means a failed `assemble` resumes without re-running the demultiplex, and each
stage gets the label and the retry that fits it. `refine` carries `process_high_memory` because
its memory is set by the number of distinct barcodes and by nothing else: a 200 GB shallow run and
a 200 GB deep one need wildly different amounts of it, and FASTQ size predicts neither.

## Wiring it in

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

## The two things that go wrong

**The barcode read is not always the payload read.** On 10x, R1 is 26 nt of cell barcode and UMI
and nothing else. `checkout` is given both mates; the later stages must then run on the mate that
carries cDNA. `payload_mate: 2` says so. Assuming R1 produces empty consensuses on every droplet
chemistry, and nothing in the run reports it as an error.

**`--rt-error` names a chemistry, not a number.** `rt` (1e-4, caps at Q40) for anything with a
reverse transcription step, `medium` (1e-5) for an ordinary polymerase and no RT, `high` (1e-6)
for a proofreading one. It is the *one-molecule* floor — 10x's Q60 requires two UMIs to agree, and
combining molecules is [arda](https://github.com/antigenomics/arda)'s job, downstream of this.

## Resources

SLURM is the pipeline's business, not the module's: these declare `label` and use `task.cpus`, and
nothing else. Every stage's output is **byte-identical at any thread count**, so a retry with
different `cpus` cannot change a result — which is what makes `errorStrategy 'retry'` with
escalating resources safe here.

```groovy
process {
    withLabel: process_high        { cpus = 16; memory = 32.GB }
    withLabel: process_medium      { cpus = 8;  memory = 16.GB }
    withLabel: process_high_memory { memory = { 64.GB * task.attempt } }
}
```

## Stubs

Every process has a `stub:` block, so `nextflow run -stub-run` exercises the whole graph without
migec installed. That is the test to run in CI before the real one.
