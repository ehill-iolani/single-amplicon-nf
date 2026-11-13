process END_TRIM {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/pandas:2.2.1"
    // with the trim on, this is the consensus that is delivered: 04_consensus/
    // gets it (MEDAKA / RACON stop publishing theirs) plus end_trim.tsv, which
    // records how many bases came off each end and what they were
    publishDir(path: { "${params.outdir}/${sample}/04_consensus" }, mode: 'copy')

    input:
    // staged under another name: the polisher's file is also called
    // <sample>.consensus.fasta, the name this step writes
    tuple val(sample), path(polished, stageAs: 'untrimmed.fasta'), path(sam)

    output:
    tuple val(sample), path("${sample}.consensus.fasta"), emit: consensus
    tuple val(sample), path("${sample}.end_trim.tsv"), emit: trim

    script:
    // sam = the consensus's own reads mapped to the polished consensus
    // (MINIMAP2_ALIGN, aliased MAP_TO_CONSENSUS); see bin/trim_consensus_ends.py
    """
    trim_consensus_ends.py \\
        --sample ${sample} \\
        --consensus ${polished} \\
        --sam ${sam} \\
        --frac ${params.end_trim_frac} \\
        --out-prefix ${sample}
    """
}
