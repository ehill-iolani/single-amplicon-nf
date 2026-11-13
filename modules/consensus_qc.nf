process CONSENSUS_QC {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    // coverage.tsv is binned depth along the consensus, small enough for the
    // frontend to plot directly
    publishDir(path: { "${params.outdir}/${sample}/05_qc" }, mode: 'copy')

    input:
    // trim = END_TRIM's end_trim.tsv, or [] with the end trim off
    tuple val(sample), path(consensus), path(paf), val(n_reads), path(trim)

    output:
    tuple val(sample), path("${sample}.qc.tsv"), emit: qc
    path "${sample}.coverage.tsv", emit: coverage

    script:
    """
    consensus_qc.py \\
        --sample ${sample} \\
        --consensus ${consensus} \\
        --paf ${paf} \\
        --n-reads ${n_reads} \\
        --min-depth ${params.min_depth} \\
        --min-primary-frac ${params.min_primary_frac} \\
        ${trim ? "--trim ${trim}" : ''} \\
        --out-prefix ${sample}
    """
}
