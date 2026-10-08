process RACON {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/racon:1.5.0--h21ec9f0_2"
    // only the final consensus is published (to 04_consensus/); with medaka on,
    // MEDAKA publishes the same filename there instead
    publishDir(path: { "${params.outdir}/${sample}/04_consensus" }, mode: 'copy', enabled: !params.enable_medaka)

    input:
    tuple val(sample), path(reads), path(draft), path(sam)

    output:
    tuple val(sample), path("${sample}.consensus.fasta"), emit: polished

    script:
    // racon always names the record "Consensus" and discards any description,
    // so the sample / read-count header is stamped on after it runs
    """
    racon -t ${task.cpus} ${reads} ${sam} ${draft} > racon.raw

    n_reads=\$(( \$(zcat ${reads} | wc -l) / 4 ))
    awk -v s="${sample}" -v n="\$n_reads" \\
        'NR == 1 { print ">" s " sample=" s " reads_used=" n " polisher=racon"; next } { print }' \\
        racon.raw > ${sample}.consensus.fasta
    """
}
