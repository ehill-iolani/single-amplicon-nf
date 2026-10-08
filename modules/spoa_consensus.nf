process SPOA_CONSENSUS {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/spoa:4.1.4--h077b44d_3"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}.draft.fasta"), emit: draft

    script:
    // one sample, one draft; the draft is only an intermediate for RACON so it
    // isn't published. The spoa container ships no seqtk -- fastq -> fasta with awk.
    """
    zcat ${reads} | awk 'NR % 4 == 1 { print ">" substr(\$0, 2) } NR % 4 == 2 { print }' > reads.fasta
    spoa reads.fasta -r 0 | awk 'NR == 1 { print ">${sample}_draft"; next } { print }' > ${sample}.draft.fasta
    """
}
