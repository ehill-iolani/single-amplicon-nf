process MEDAKA {
    tag "$sample"
    label 'process_medium'
    // ONT's own image, pinned to the same digest as ulana-nf: medaka 2.2.1,
    // which knows the v5.0.0 models (the v1.11.3 tag does not). Not biocontainers:
    // that build is amd64-only and SIGILLs under Docker emulation on Apple Silicon.
    container "${params.dockerhub_registry}/ontresearch/medaka:shaf39439188c323053e66cf79618fae1ab33d1b38a"
    publishDir(path: { "${params.outdir}/${sample}/04_consensus" }, mode: 'copy')

    input:
    tuple val(sample), path(reads), path(racon_fasta)

    output:
    tuple val(sample), path("${sample}.consensus.fasta"), emit: consensus

    script:
    """
    medaka_consensus -i ${reads} -d ${racon_fasta} -o medaka_out -t ${task.cpus} -m ${params.medaka_model}

    n_reads=\$(( \$(zcat ${reads} | wc -l) / 4 ))
    awk -v s="${sample}" -v n="\$n_reads" \\
        'NR == 1 { print ">" s " sample=" s " reads_used=" n " polisher=medaka"; next } { print }' \\
        medaka_out/consensus.fasta > ${sample}.consensus.fasta
    """
}
