process CUTADAPT {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/cutadapt:4.9--py310h1fe012e_3"
    publishDir(path: { "${params.outdir}/${sample}/02_trimmed" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.trimmed.fastq.gz"), emit: trimmed

    script:
    // Primers are given 5'->3' as ordered. They are joined into one linked
    // adapter, FWD...revcomp(REV), and --revcomp makes cutadapt flip reads that
    // match on the reverse strand, so every trimmed read comes out in the
    // forward-primer orientation (which spoa needs). Reads where the primers
    // are not found are discarded -- they are off-target by definition.
    // Without primers the reads pass through untouched, and orientation is
    // left to the purity check (DOMINANT_CLUSTER).
    if (params.fwd_primer && params.rev_primer)
        """
        fwd=\$(echo "${params.fwd_primer}" | tr 'a-z' 'A-Z')
        rev_rc=\$(echo "${params.rev_primer}" | tr 'a-z' 'A-Z' | rev | tr 'ACGTUMRWSYKVHDBN' 'TGCAAKYWSRMBDHVN')
        cutadapt -j ${task.cpus} --revcomp --discard-untrimmed \\
            -g "\${fwd}...\${rev_rc}" \\
            -o ${sample}.trimmed.fastq.gz ${fastq}
        """
    else
        """
        cp ${fastq} ${sample}.trimmed.fastq.gz
        """
}
