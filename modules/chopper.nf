process CHOPPER {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/chopper:0.7.0--hdcf5f25_0"
    publishDir(path: { "${params.outdir}/${sample}/01_filtered" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.filtered.fastq.gz"), emit: filtered

    script:
    // max_len is optional ('' from the platform's form counts as unset)
    def maxlen = params.max_len ? "--maxlength ${params.max_len}" : ""
    """
    set -o pipefail

    zcat -f ${fastq} \\
      | chopper -t ${task.cpus} -q ${params.min_qual} -l ${params.min_len} ${maxlen} \\
      | gzip > ${sample}.filtered.fastq.gz
    """
}
