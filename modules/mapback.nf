process MAPBACK {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/minimap2:2.28--he4a0461_3"

    input:
    tuple val(sample), path(consensus), path(reads)

    output:
    tuple val(sample), path(consensus), path("${sample}.paf"), env('N_READS'), emit: paf

    script:
    // base-level (-c) PAF so the match count is exact. Per-position depth and
    // identity are worked out from the PAF in CONSENSUS_QC (stdlib python), so
    // neither samtools nor mosdepth is needed.
    """
    minimap2 -c -x map-ont -t ${task.cpus} ${consensus} ${reads} > ${sample}.paf
    N_READS=\$(( \$(zcat ${reads} | wc -l) / 4 ))
    """
}
