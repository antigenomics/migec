# migec as nf-core-style local modules

```
main.nf                                    a runnable entry point: --mode consensus | ctdna | airr
nextflow.config                            defaults, every one read with getOrDefault

modules/migec/checkout/main.nf             reads     -> tagged FASTQ + QC tables
modules/migec/refine/main.nf               tagged    -> corrected FASTQ + barcode table + cell calls
modules/migec/assemble/main.nf             corrected -> one consensus per molecule
modules/migec/plot/main.nf                 the tables -> SVG figures (needs gnuplot)

modules/downstream/align/main.nf           consensus -> BAM, tags carried and CHECKED
modules/downstream/callvariants/main.nf    BAM       -> VCF (LoFreq or Mutect2)
modules/downstream/arda/main.nf            consensus -> AIRR clonotypes

subworkflows/migec/main.nf                 the three migec stages chained
subworkflows/migec_ctdna/main.nf           + align + call: rare somatic variants
subworkflows/migec_airr/main.nf            + arda: immune repertoires
```

## Run it

```bash
# rare variants from a ctDNA panel
nextflow run integrations/nextflow -profile docker \
    --mode ctdna --input 'data/*_R{1,2}.fq.gz' --preset tso500 \
    --fasta ref.fa --outdir results/

# an immune repertoire
nextflow run integrations/nextflow --mode airr --input 'data/*_R{1,2}.fq.gz' --preset migec

# just the consensus, to hand to something else
nextflow run integrations/nextflow --mode consensus --input 'data/*_R{1,2}.fq.gz' --preset 10x-v2

# walk the whole graph with nothing installed
nextflow run integrations/nextflow --mode ctdna --input 'x_R{1,2}.fq.gz' --fasta r.fa -stub-run
```

## Or include the subworkflows in your own pipeline

```groovy
include { MIGEC_CTDNA } from './integrations/nextflow/subworkflows/migec_ctdna/main'

workflow {
    ch_reads = Channel.fromFilePairs(params.input)
        .map { id, files -> [ [ id: id, preset: 'tso500', payload_mate: 1 ], files ] }

    MIGEC_CTDNA(ch_reads, Channel.value(file(params.fasta)),
                Channel.value(file("${params.fasta}.fai")), file('NO_FILE'))

    MIGEC_CTDNA.out.vcf.view()
}
```

Per-sample keys in `meta` beat the `params.*` defaults, so one run can mix chemistries:
`bc_pattern`, `preset`, `read_structure`, `read_structure2`, `barcodes`, `max_offset`,
`payload_mate`, `expect_cells`, `rt_error`, `contig`, `fast`, `aligner`, `caller`, `species`.

## The two rules the downstream modules encode

**Collapse first, then align once.** Aligning raw reads and grouping on *(position, UMI)* is the
other order in use — fgbio, UMI-tools, UMIErrorCorrect — and it costs one alignment per *read*
rather than per *molecule*, with the aligner seeing uncorrected sequence.

**A standard variant caller, never a UMI-aware one.** LoFreq, Mutect2, FreeBayes and VarDict read
a BAM and ignore `RX`, so after `assemble` their depth already is a molecule count. UMI-VarCal and
UMIErrorCorrect group and consense themselves, so they *replace* `assemble` rather than following
it. Never set a family-size filter downstream of `assemble`: every family has size 1 by
construction, so `--min-family-size 3` discards the whole library and reports zero variants without
an error.

`docs/variants.rst` has the full table and the measurement behind the LoFreq default.

## What is verified and what is not

The align module **checks its own output**: if no `MI:Z:` tag survives into the BAM it exits
non-zero, naming the flag each aligner needs (`-y` for minimap2 and `minibwa map`, `-C` for
`bwa mem`). That check exists because the failure is otherwise silent and only shows up as an
inexplicably untagged BAM much later.

Never: **nextflow is not installed on the machine these modules were written on**, so they are
reviewed against the nf-core module spec rather than verified by a run. Treat your first
`-stub-run` as the verification — every process has a `stub:` block for exactly that. The migec
commands inside them *are* tested, by `tests/` and by `integrations/slurm/`, which runs the same
three stages as ordinary shell.

**The full write-up is [`docs/nextflow.rst`](../../docs/nextflow.rst)** — why it is three processes
and not one, the resource labels, and the three things that go wrong (the payload mate, a
per-sample `false`, and `--rt-error` naming a chemistry rather than a number).
