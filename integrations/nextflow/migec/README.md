# migec as an nf-core-style module

Drops into any pipeline that hands you FASTQ pairs — [nf-core/airrflow](https://nf-co.re/airrflow)
being the obvious one, since it does everything downstream of a consensus and nothing upstream of
it.

```
modules/local/migec/
├── main.nf
├── meta.yml
├── nextflow.config
└── environment.yml
```

```groovy
include { MIGEC } from '../modules/local/migec/main'

workflow {
    MIGEC(ch_reads)          // tuple val(meta), path(reads)
    MIGEC.out.consensus.view()
}
```

```groovy
// nextflow.config
includeConfig 'modules/local/migec/nextflow.config'
params.migec_bc_pattern = 'XXXXXXXXXXXXXXXXNNNNNNNNNN'   // 10x 5' v2
params.migec_max_offset = 0
params.migec_payload_mate = 2
```

## The two settings that are not conveniences

**`migec_bc_pattern` has no default.** The layout is a property of the library prep, and a wrong
one produces barcodes rather than an error. If it is genuinely unknown, `migec suggest reads.fq.gz`
reads it off the data — it recovered a 9 nt Primer ID and its 23 nt anchor from a public HIV run
with nothing supplied but the FASTQ.

**`migec_payload_mate` is 2 on every droplet chemistry.** The 10x barcode read is 26 nt of cell
barcode and UMI and *nothing else*; the cDNA is all on the other mate. Run `refine`/`assemble` on
R1 there and you get empty consensuses — migec reports a payload clonality of 1.0 saying so, which
is easy to miss in a pipeline log and impossible to miss here.

## Resources

Memory scales with **distinct barcodes**, never with reads:

| stage | holds | threads |
|---|---|---|
| `checkout` | chunk × threads, plus the UMI counters | `task.cpus`, output byte-identical at any count |
| `refine` | the barcode table, ~96 B per distinct barcode; reads streamed three times | 1 |
| `assemble` | one range-partitioned bucket | 1 |

Measured: 3.16 M 10x VDJ reads over 311 k barcodes peaked **658 MB**. Budget by expected molecules,
not by FASTQ size. `label 'process_medium'` suits a few million reads; raise it for a run with
10⁸ distinct barcodes, where the `refine` table alone is ~10 GB.

## SLURM

Nothing migec-specific — the module declares `label`, `task.cpus` and nothing else, so the
executor is entirely your pipeline's business:

```groovy
process {
    executor = 'slurm'
    queue    = 'normal'
    withName: 'MIGEC' {
        cpus   = 16
        memory = { 8.GB * task.attempt }
        time   = { 2.h * task.attempt }
    }
}
```

⚠ `checkout` is the only stage that threads; `refine` and `assemble` are single-threaded by
construction, so a 64-core allocation buys nothing after the first stage. Ask for the cores
`checkout` can use and no more.

⚠ **On a shared cluster, never run this on the login node.** Submit it.

## Without Nextflow

The three commands are the whole interface, and there is nothing in the module a shell loop cannot
do:

```bash
migec checkout R1.fq.gz R2.fq.gz --bc-pattern "$PATTERN" --max-offset 0 -t "$SLURM_CPUS_ON_NODE" -o co/
migec refine    co/sample_R2.fq.gz -o ref/
migec assemble  ref/sample.fq.gz   -o asm/
```
