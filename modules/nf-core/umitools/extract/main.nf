process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_long"

    conda (params.enable_conda ? "bioconda::umi_tools=1.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/umi_tools:1.1.2--py38h4a8c8d9_0' :
        'quay.io/biocontainers/umi_tools:1.1.2--py38h4a8c8d9_0' }"

    input:
    tuple val(meta), path(reads)
    path whitelist

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    umi_tools \\
        extract \\
        --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \\
        -I ${reads[0]} \\
        --read2-in=${reads[1]} \\
        -S ${prefix}.umi_extract_1.fastq.gz \\
        --read2-out=${prefix}.umi_extract_2.fastq.gz \\
        --whitelist=$whitelist \\
        $args \\
        > ${prefix}.umi_extract.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
    END_VERSIONS
    """
}
