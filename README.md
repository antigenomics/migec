> Deprecated. This is the legacy version, kept on the [`legacy-v1`](https://github.com/antigenomics/migec/tree/legacy-v1) branch and no longer supported; use the latest version at [antigenomics/migec](https://github.com/antigenomics/migec).

[![JitPack](https://jitpack.io/v/mikessh/migec.svg)](https://jitpack.io/#mikessh/migec)

## MiGEC: Molecular Identifier Guided Error Correction pipeline  

This pipeline provides several useful tools for analysis of immune repertoire sequencing data. Its main feature is the ability to use information from unique nucleotide tags (UMIs, see this [paper](http://www.nature.com/nmeth/journal/v9/n1/full/nmeth.1778.html) for details), which are attached to molecules before sequencing library preparation and allow to backtrack the original sequence of molecule. UMIs make it possible to computationally filter nearly all experimental errors from resulting immune receptor sequences.

This pipeline was designed for libraries sequenced using Illumina MiSeq and HiSeq and the main requirement for sequencing reads is that they should contain the entire CDR3 region of immune receptor gene. Sequencing libraries with high over-sequencing, i.e. ones that have 5+ reads per starting molecule (unique UMI tag), should be used for optimal error elimination.

Several modules of the pipeline, such as de-multiplexing and CDR3 extraction could be utilized for a wider range of datasets.

Compiled binaries are available from [here](https://github.com/mikessh/migec/releases/latest). You can download them and execute as

```bash
java -jar migec.jar ...
```

Make sure that you've specified the full/correct path to jar file. In case of Java Heap Space exception, you can increase the JVM memory limit by adding ``-Xmx20G`` (for extra 20G) after the ``-jar`` argument.

The software is cross-platform and requires Java v1.8+ to run.

Easy installation on MacOS/Linux via [Homebrew](http://brew.sh/) or [Linuxbrew](http://linuxbrew.sh/):
```bash
brew tap homebrew/science
brew tap mikessh/repseq
brew install migec
migec Checkout ...
```
See [homebrew-repseq](https://github.com/mikessh/homebrew-repseq) for other RepSeq analysis software Homebrew installers.

For more details please see the [paper](http://www.nature.com/nmeth/journal/v11/n6/abs/nmeth.2960.html) describing MiGEC.

Full documentation is provided via [ReadTheDocs](http://migec.readthedocs.org/en/latest/index.html). You might be also interested in taking the following [tutorial](http://repseq-tutorial.readthedocs.org/en/latest/).

Please cite the tool as:

> Shugay M *et al.* Towards error-free profiling of immune repertoires. Nature Methods **11**, 653–655 (2014)

You may be also interested in the following [Nature Protocol](https://www.nature.com/articles/nprot.2016.093) describing MIGEC applications for full-length immunoglobulin sequencing.

### Some notes on MiGEC study

The basic idea behind this project is to incorporate unique molecular identifier tags to trace individual cDNA molecules during library preparation and sequencing in order to eliminate experimental errors. This is especially crucial for B-cell studies as it is necessary to tell real hypermutations from hot-spot PCR and sequencing errors.

Here we provide a modified protocol for performing nearly error-free high-throughput sequencing of TCR and antibody repertoires, as well as a novel error-correction algorithm that successfully handles hot-spot PCR errors.

A detailed description of the study and protocol used could be found in our [Nature Methods paper](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.2960.html).

A software pipeline that allows analysis of immune repertoire sequencing data (for both T- and B-cells) prepared using our protocol is available at GitHub. Here are the links to readme and a package containing the executable Jar. The pipeline is written using Java/Groovy and has moderate system requirements for large datasets: it will analyze a HiSeq lane in 4-5 hours on a commodity server with 8 core Intel Xeon and 36 Gb RAM.

Raw datasets used in the study are deposed in SRA under accessions [SRP040329](http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP040329). The barcode sequences that are necessary for de-multiplexing and UMI extraction with Checkout utility are available [here](https://github.com/mikessh/migec/blob/master/misc/barcodes.txt?raw=true). RNA for immunoglobulin heavy chain (IGH) from more than 600,000 B-cells from blood of a healthy individual was sequenced, yielding a repertoire of more than 100,000 clonotypes.

A set of spike-ins EHEB, EHEB-V1 (with one mismatch difference from EHEB-V0) and EHEB-V2 (with two mismatches) were added in order to monitor the efficiency of error filtration.

Have a glance at the analyzed data (done with legacy MIGEC version):

* The spreadsheet demonstrating efficient error elimination for spike-in clonotypes with known sequences [[Download](https://github.com/mikessh/migec/blob/master/misc/Exp2-spikein-table.xlsx?raw=true)]
* The spreadsheet with whole repertoire [[Download](https://github.com/mikessh/migec/blob/master/misc/Exp2-all-cdr-migec.xlsx?raw=true)]
* An interactive network of B-cell clonal trees (created using [Cytoscape](http://cytoscape.org/)) [[Download](https://github.com/mikessh/migec/blob/master/misc/Exp2-all-cdr-migec.cys?raw=true)]
