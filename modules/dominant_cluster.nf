process DOMINANT_CLUSTER {
    tag "$sample"
    label 'process_medium'
    container "${params.container_registry}/biocontainers/vsearch:2.30.6--h0bb26bb_0"
    publishDir(path: { "${params.outdir}/${sample}/03_cluster" }, mode: 'copy')

    input:
    tuple val(sample), path(fastq)

    output:
    tuple val(sample), path("${sample}.dominant.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.purity.tsv"), emit: purity

    script:
    /*
     * Purity check for a sample that is supposed to hold ONE sequence: cluster
     * the reads, keep the largest cluster, and report what share of the reads
     * it holds. A clean sample lands ~all reads in one cluster; a mixed or
     * contaminated one doesn't, and the consensus is then built from the
     * majority sequence instead of from a blend of several (which spoa would
     * turn into a chimera).
     *
     * cluster_id is deliberately loose (0.70): it only has to tolerate raw ONT
     * error between two reads of the SAME sequence, not separate species. Two
     * reads that are each ~6% wrong differ from each other by ~12% (vsearch
     * counts indels as mismatches), so reads sit at ~80-88% identity to their
     * centroid and 0.80 splits a clean sample into several clusters (measured
     * on the synthetic test data: clean sample 1.00 dominant at 0.75 and below,
     * 0.74 at 0.80). Retune per basecaller/chemistry. The flip side: this
     * catches off-target/unrelated contamination, NOT closely related
     * co-amplified variants (<~15% divergent) -- those fall in one cluster and
     * show up as low identity / ambiguous bases in CONSENSUS_QC instead.
     *
     * --strand both: reads come off the pore in either orientation. Members that
     * hit the centroid on the minus strand are reverse-complemented (qualities
     * reversed) so the whole cluster shares one orientation.
     */
    """
    set -o pipefail

    zcat ${fastq} > input.fastq

    # fastq -> fasta for clustering (vsearch's fastq parser would also reject
    # modern Q-scores above its default --fastq_qmax of 41)
    awk 'NR % 4 == 1 { print ">" substr(\$0, 2) } NR % 4 == 2 { print }' input.fastq > reads.fasta

    vsearch --cluster_fast reads.fasta \\
        --id ${params.cluster_id} \\
        --strand both \\
        --threads ${task.cpus} \\
        --uc clusters.uc

    # clusters.uc: S = centroid, H = member hit; field 2 = cluster number,
    # field 5 = strand vs. the centroid, field 9 = read id.
    # Dominant cluster = most members; ties go to the lower cluster number so
    # the choice never depends on awk's hash order.
    awk -F'\\t' '
        \$1 == "S" || \$1 == "H" { n[\$2]++; total++ }
        END {
            best = -1; bn = 0; second = 0; k = 0
            for (c in n) {
                k++
                if (n[c] > bn || (n[c] == bn && c + 0 < best + 0)) {
                    if (bn > second) second = bn
                    bn = n[c]; best = c
                } else if (n[c] > second) {
                    second = n[c]
                }
            }
            printf "%s\\t%d\\t%d\\t%d\\t%d\\n", best, total, k, bn, second
        }
    ' clusters.uc > dominant.txt
    read dom total k dom_n second_n < dominant.txt

    printf 'sample\\treads_clustered\\tn_clusters\\tdominant_reads\\tdominant_fraction\\tsecond_fraction\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample}" "\$total" "\$k" "\$dom_n" \\
        "\$(awk -v a=\$dom_n -v b=\$total 'BEGIN { printf "%.4f", b ? a / b : 0 }')" \\
        "\$(awk -v a=\$second_n -v b=\$total 'BEGIN { printf "%.4f", b ? a / b : 0 }')" \\
        > ${sample}.purity.tsv

    # pull the dominant cluster's reads out of the original fastq (qualities
    # intact, racon/medaka need them), split by strand for the flip below
    awk -F'\\t' -v dom=\$dom '
        NR == FNR { if ((\$1 == "S" || \$1 == "H") && \$2 == dom) st[\$9] = \$5; next }
        FNR % 4 == 1 { id = substr(\$0, 2) }
        FNR % 4 == 2 { seq = \$0 }
        FNR % 4 == 0 {
            if (id in st) {
                out = (st[id] == "-") ? "minus.tsv" : "plus.tsv"
                print id "\\t" seq "\\t" \$0 > out
            }
        }
    ' clusters.uc input.fastq
    touch plus.tsv minus.tsv

    cut -f1 minus.tsv > minus.id
    cut -f2 minus.tsv | rev | tr 'ACGTacgt' 'TGCAtgca' > minus.seq
    cut -f3 minus.tsv | rev > minus.qual
    paste minus.id minus.seq minus.qual > minus.oriented.tsv

    cat plus.tsv minus.oriented.tsv \\
        | awk -F'\\t' '{ printf "@%s\\n%s\\n+\\n%s\\n", \$1, \$2, \$3 }' \\
        | gzip > ${sample}.dominant.fastq.gz

    echo "DOMINANT_CLUSTER ${sample}: \$dom_n of \$total reads in the dominant cluster (\$k clusters at --id ${params.cluster_id})" >&2
    """
}
