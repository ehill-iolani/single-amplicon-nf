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
    // only when the second-largest cluster is big enough to build a consensus
    // from (see below); otherwise nothing is emitted for the sample
    tuple val(sample), path("${sample}.secondary.fastq.gz"), optional: true, emit: secondary

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
     *
     * Second sequence: when the second-largest cluster holds at least
     * --min_secondary_frac of the clustered reads (and at least --min_reads
     * reads, the same floor a sample needs for a consensus at all), its reads
     * are written out the same way, as <sample>.secondary.fastq.gz, so the
     * workflow can build a consensus from them too. That is the "two real
     * sequences in one barcode" case; a scatter of small clusters is noise and
     * gets nothing. Off with --enable_secondary_consensus false.
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
    # the choice never depends on awk's hash order. Second = the biggest of the
    # rest, same tie rule (-1 and 0 reads when there is only one cluster).
    awk -F'\\t' '
        \$1 == "S" || \$1 == "H" { n[\$2]++; total++ }
        END {
            best = -1; bn = 0; next_id = -1; sn = 0; k = 0
            for (c in n) {
                k++
                if (n[c] > bn || (n[c] == bn && c + 0 < best + 0)) { best = c; bn = n[c] }
            }
            for (c in n) {
                if (c == best) continue
                if (n[c] > sn || (n[c] == sn && c + 0 < next_id + 0)) { next_id = c; sn = n[c] }
            }
            printf "%s\\t%d\\t%d\\t%d\\t%d\\t%s\\n", best, total, k, bn, sn, next_id
        }
    ' clusters.uc > dominant.txt
    read dom total k dom_n second_n second_id < dominant.txt

    printf 'sample\\treads_clustered\\tn_clusters\\tdominant_reads\\tdominant_fraction\\tsecond_fraction\\n%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample}" "\$total" "\$k" "\$dom_n" \\
        "\$(awk -v a=\$dom_n -v b=\$total 'BEGIN { printf "%.4f", b ? a / b : 0 }')" \\
        "\$(awk -v a=\$second_n -v b=\$total 'BEGIN { printf "%.4f", b ? a / b : 0 }')" \\
        > ${sample}.purity.tsv

    # pull one cluster's reads out of the original fastq (qualities intact,
    # racon/medaka need them), split by strand for the flip below:
    #   extract_cluster <cluster number> <output fastq.gz>
    extract_cluster() {
        local cluster="\$1" out="\$2"
        awk -F'\\t' -v dom="\$cluster" '
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
            | gzip > "\$out"
        rm -f plus.tsv minus.tsv minus.id minus.seq minus.qual minus.oriented.tsv
    }

    extract_cluster "\$dom" ${sample}.dominant.fastq.gz

    if [ "${params.enable_secondary_consensus}" = "true" ] && awk -v n=\$second_n -v t=\$total -v min_n=${params.min_reads} -v min_f=${params.min_secondary_frac} \\
            'BEGIN { exit !(t > 0 && n >= min_n && n / t >= min_f) }'; then
        extract_cluster "\$second_id" ${sample}.secondary.fastq.gz
        echo "DOMINANT_CLUSTER ${sample}: second cluster has \$second_n of \$total reads; building a second consensus from it" >&2
    fi

    echo "DOMINANT_CLUSTER ${sample}: \$dom_n of \$total reads in the dominant cluster (\$k clusters at --id ${params.cluster_id})" >&2
    """
}
