process CONSENSUS_QC {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    // coverage.tsv is binned depth along the consensus, small enough for the
    // frontend to plot directly
    publishDir(path: { "${params.outdir}/${sample}/05_qc" }, mode: 'copy')

    input:
    tuple val(sample), path(consensus), path(paf), val(n_reads)

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
        --out-prefix ${sample}
    """
}
