process DOWNSAMPLE {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/seqkit:2.8.2--h9ee0642_1"

    input:
    tuple val(sample), path(fastq)

    output:
    // n_used is read by the workflow to decide whether the sample is usable
    tuple val(sample), path("${sample}.reads.fastq.gz"), env('N_USED'), emit: reads
    tuple val(sample), path("${sample}.counts.tsv"), emit: counts

    script:
    // Read headers are cut to the bare id first: Dorado headers carry
    // tab-separated tags (qs:f:.. mx:i:..) that racon/vsearch/minimap2 each
    // treat a little differently, and nothing downstream needs them.
    //
    // Reads are then sorted by id so the rest of the pipeline is reproducible.
    // CHOPPER is multithreaded and emits the same reads in a different order each
    // run, which would make the seeded sample below pick different reads, and
    // spoa/vsearch results depend on read order too. Each 4-line read is joined
    // onto one line (paste), sorted by id (byte order, so the locale can't change
    // it), and split again. `sample -2 -n` is seqkit's exact two-pass sampler,
    // which needs a file, not a pipe.
    """
    set -o pipefail

    seqkit seq -i ${fastq} -o ids.fastq.gz
    zcat ids.fastq.gz | paste - - - - | LC_ALL=C sort -k1,1 | tr '\\t' '\\n' | gzip > sorted.fastq.gz
    N_IN=\$(( \$(zcat sorted.fastq.gz | wc -l) / 4 ))

    if [ "${params.downsample_reads}" -gt 0 ] && [ "\$N_IN" -gt "${params.downsample_reads}" ]; then
        seqkit sample -2 -n ${params.downsample_reads} -s ${params.seed} sorted.fastq.gz -o ${sample}.reads.fastq.gz
    else
        mv sorted.fastq.gz ${sample}.reads.fastq.gz
    fi
    N_USED=\$(( \$(zcat ${sample}.reads.fastq.gz | wc -l) / 4 ))

    printf 'sample\\treads_after_trim\\treads_used\\n%s\\t%s\\t%s\\n' "${sample}" "\$N_IN" "\$N_USED" > ${sample}.counts.tsv
    """
}
