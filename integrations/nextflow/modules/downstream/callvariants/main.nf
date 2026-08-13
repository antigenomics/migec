// Call rare variants from the consensus BAM.
//
// Never: USE A STANDARD CALLER HERE, NOT A UMI-AWARE ONE. LoFreq, Mutect2, FreeBayes and VarDict
// read a BAM and never look at `RX` -- so after `assemble` their depth already IS a molecule count
// and their allele fraction already IS a molecule fraction. UMI-VarCal and UMIErrorCorrect do
// their own grouping and consensus, which means they REPLACE assemble rather than following it;
// running both counts every molecule once and then collapses the result again.
//
// Never: DO NOT SET A FAMILY-SIZE FILTER. After assemble every family has size 1 by construction,
// so a `--min-family-size 3` discards the entire library and reports zero variants without an
// error. docs/variants.rst has the table of which caller sits on which side of this line.
//
// The default is LoFreq: on the only independent ctDNA comparison of the six (Maruzani et al.
// 2024, doi:10.1186/s12864-024-10737-w) it returned the fewest false positives of any standard
// caller at every depth, and was second only to Mutect2 on sensitivity.

process DOWNSTREAM_CALLVARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.getOrDefault('caller_container',
                                  'quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4')

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.vcf.gz"), emit: vcf
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def caller = meta.caller ?: params.getOrDefault('downstream_caller', 'lofreq')
    def args   = params.getOrDefault('downstream_caller_args', '')
    if (caller == 'lofreq') {
        """
        lofreq call-parallel --pp-threads ${task.cpus} \\
            -f ${fasta} -o ${prefix}.vcf ${args} ${bam}
        bgzip -f ${prefix}.vcf && tabix -p vcf ${prefix}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            lofreq: \$(lofreq version 2>&1 | head -1 | awk '{print \$2}')
        END_VERSIONS
        """
    }
    else if (caller == 'mutect2') {
        // FilterMutectCalls is not optional. Mutect2's own authors require it, and the published
        // comparison that scored Mutect2 as low-specificity deliberately omitted it.
        """
        gatk Mutect2 -R ${fasta} -I ${bam} -O ${prefix}.unfiltered.vcf.gz ${args}
        gatk FilterMutectCalls -R ${fasta} \\
            -V ${prefix}.unfiltered.vcf.gz -O ${prefix}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk: \$(gatk --version 2>&1 | grep -o 'v[0-9.]*' | head -1)
        END_VERSIONS
        """
    }
    else {
        error "downstream_caller must be lofreq or mutect2; got '${caller}'. " +
              "UMI-aware callers replace assemble rather than following it -- see docs/variants.rst"
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        caller: stub
    END_VERSIONS
    """
}
