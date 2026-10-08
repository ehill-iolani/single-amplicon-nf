process DOWNSAMPLE {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/seqkit:2.8.2--h9ee0642_1"

    input:
    tuple val(sample), path(fastq)

    output:
    // n_used is read by the workflow to decide whether the sample is usable
    tuple val(sample), path("${sample}.reads.fastq.gz"), env(N_USED), emit: reads
    tuple val(sample), path("${sample}.counts.tsv"), emit: counts

    script:
    // Read headers are cut to the bare id first: Dorado headers carry
    // tab-separated tags (qs:f:.. mx:i:..) that racon/vsearch/minimap2 each
    // treat a little differently, and nothing downstream needs them.
    // `sample -2 -n` is seqkit's exact two-pass mode (needs a file, not a
    // pipe); the fixed --seed makes a rerun reproduce the same consensus.
    """
    set -o pipefail

    seqkit seq -i ${fastq} -o ids.fastq.gz
    N_IN=\$(( \$(zcat ids.fastq.gz | wc -l) / 4 ))

    if [ "${params.downsample_reads}" -gt 0 ] && [ "\$N_IN" -gt "${params.downsample_reads}" ]; then
        seqkit sample -2 -n ${params.downsample_reads} -s ${params.seed} ids.fastq.gz -o ${sample}.reads.fastq.gz
    else
        mv ids.fastq.gz ${sample}.reads.fastq.gz
    fi
    N_USED=\$(( \$(zcat ${sample}.reads.fastq.gz | wc -l) / 4 ))

    printf 'sample\\treads_after_trim\\treads_used\\n%s\\t%s\\t%s\\n' "${sample}" "\$N_IN" "\$N_USED" > ${sample}.counts.tsv
    """
}
