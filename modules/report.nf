process BUILD_REPORT {
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    publishDir "${params.outdir}/final_report", mode: 'copy'

    input:
    path count_files
    path purity_files
    path qc_files
    path consensus_fastas

    output:
    path "consensus_summary.tsv", emit: summary
    path "all_consensus.fasta"
    path "run_qc_summary.html"
    // only when a sample got a second consensus
    path "secondary_consensus.tsv", optional: true
    path "secondary_consensus.fasta", optional: true

    script:
    // purity_files, qc_files and consensus_fastas can each be empty (check off /
    // no sample got a consensus); every
    // multi-value option below is followed by another flag so an empty list
    // stays a valid command line
    """
    build_report.py \\
        --counts ${count_files} \\
        --purity ${purity_files} \\
        --qc ${qc_files} \\
        --consensus ${consensus_fastas} \\
        --min-reads ${params.min_reads} \\
        --min-dominant-frac ${params.min_dominant_frac} \\
        --secondary-suffix ${params.secondary_suffix} \\
        --out-summary consensus_summary.tsv \\
        --out-fasta all_consensus.fasta \\
        --out-html run_qc_summary.html
    """
}
