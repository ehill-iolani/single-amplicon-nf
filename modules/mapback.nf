process MAPBACK {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/minimap2:2.28--he4a0461_3"
    // Only the alignments are published: the pattern keeps the consensus and the
    // PAF (both also outputs, for CONSENSUS_QC) out of 05_qc/, which would
    // otherwise get every output.
    publishDir(path: { "${params.outdir}/${sample}/05_qc" }, mode: 'copy', pattern: '*.alignments.sam.gz')

    input:
    tuple val(sample), path(consensus), path(reads)

    output:
    tuple val(sample), path(consensus), path("${sample}.paf"), env('N_READS'), emit: paf
    tuple val(sample), path("${sample}.alignments.sam.gz"), emit: sam

    script:
    // base-level (-c) PAF so the match count is exact. Per-position depth and
    // identity are worked out from the PAF in CONSENSUS_QC (stdlib python), so
    // neither samtools nor mosdepth is needed.
    //
    // The same mapping is also written as SAM, for the frontend's read viewer
    // (it draws each read against the consensus and marks the bases that
    // differ). Same minimap2 settings, so it shows the alignments the QC
    // numbers came from. It is made for display, so it is trimmed to stay
    // small: unmapped reads are dropped, base qualities are blanked, and the
    // @PG line (the command line, with work-dir paths) is left out. It is not
    // sorted or indexed. --secondary=no keeps one line per alignment.
    // (The awk is portable -- no gawk-only and() -- since the container's awk
    // is not guaranteed to be gawk. FLAG 4 = unmapped.)
    """
    set -o pipefail

    minimap2 -c -x map-ont -t ${task.cpus} ${consensus} ${reads} > ${sample}.paf
    N_READS=\$(( \$(zcat ${reads} | wc -l) / 4 ))

    minimap2 -a -x map-ont --secondary=no -t ${task.cpus} ${consensus} ${reads} \\
        | awk 'BEGIN { OFS = "\\t" }
               /^@PG/ { next }
               /^@/   { print; next }
               int(\$2 / 4) % 2 == 1 { next }
               { \$11 = "*"; print }' \\
        | gzip -c > ${sample}.alignments.sam.gz
    """
}
