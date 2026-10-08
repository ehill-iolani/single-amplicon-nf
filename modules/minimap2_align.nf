process MINIMAP2_ALIGN {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/minimap2:2.28--he4a0461_3"

    input:
    tuple val(sample), path(reads), path(draft)

    output:
    tuple val(sample), path(reads), path(draft), path("${sample}.sam"), emit: aligned

    script:
    // the SAM is only an intermediate for RACON, so it isn't published
    """
    minimap2 -ax map-ont -t ${task.cpus} ${draft} ${reads} > ${sample}.sam
    """
}
