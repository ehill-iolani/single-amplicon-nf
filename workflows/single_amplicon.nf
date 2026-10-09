include { MERGE_FASTQ      } from '../modules/merge_fastq.nf'
include { CHOPPER           } from '../modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from '../modules/read_stats.nf'
include { CUTADAPT          } from '../modules/cutadapt.nf'
include { DOWNSAMPLE        } from '../modules/downsample.nf'
include { DOMINANT_CLUSTER  } from '../modules/dominant_cluster.nf'
include { SPOA_CONSENSUS    } from '../modules/spoa_consensus.nf'
include { MINIMAP2_ALIGN    } from '../modules/minimap2_align.nf'
include { RACON             } from '../modules/racon.nf'
include { MEDAKA            } from '../modules/medaka.nf'
include { MAPBACK           } from '../modules/mapback.nf'
include { CONSENSUS_QC      } from '../modules/consensus_qc.nf'
include { BUILD_REPORT      } from '../modules/report.nf'

workflow SINGLE_AMPLICON {

    take:
    reads_ch    // tuple(sample, fastq)

    main:
    // 1. merge multi-part fastq(.gz) per barcode, then length/quality filter
    MERGE_FASTQ(reads_ch)
    CHOPPER(MERGE_FASTQ.out.merged)

    // read length / Q-score summary before vs. after filtering -- a side
    // branch, nothing downstream depends on it
    if (params.enable_read_stats) {
        READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
        READ_STATS_REPORT(
            READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
            READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
        )
    }

    // 2. primer trimming (+ orientation, when primers are given) and
    //    random downsampling to --downsample_reads
    CUTADAPT(CHOPPER.out.filtered)
    DOWNSAMPLE(CUTADAPT.out.trimmed)

    // 3. a sample with too few reads gets no consensus but is still reported
    //    (FAIL, too_few_reads) -- it must not stop the other samples
    usable_ch = DOWNSAMPLE.out.reads
        .filter { sample, fastq, n_used -> n_used.toInteger() >= params.min_reads }
        .map { sample, fastq, n_used -> tuple(sample, fastq) }
    DOWNSAMPLE.out.reads
        .filter { sample, fastq, n_used -> n_used.toInteger() < params.min_reads }
        .subscribe { sample, fastq, n_used ->
            log.warn "Sample '${sample}': ${n_used} reads after filtering is below --min_reads ${params.min_reads}; no consensus will be built"
        }

    // 4. purity check: keep only the dominant read cluster (strand-aware, so
    //    every kept read also shares the centroid's orientation, which spoa needs).
    //    When the second-largest cluster is big enough to be a second sequence
    //    (--min_secondary_frac), its reads come out too, under the id
    //    "<sample><secondary_suffix>"; they take the same steps below as a sample
    //    of their own, so the second consensus gets the same polish and QC
    if (params.enable_purity_check) {
        DOMINANT_CLUSTER(usable_ch)
        consensus_reads = DOMINANT_CLUSTER.out.reads
        purity_files    = DOMINANT_CLUSTER.out.purity.map { sample, tsv -> tsv }.collect().ifEmpty([])
        secondary_reads = DOMINANT_CLUSTER.out.secondary.map { sample, fastq -> tuple("${sample}${params.secondary_suffix}".toString(), fastq) }
    } else {
        consensus_reads = usable_ch
        purity_files    = Channel.value([])
        secondary_reads = Channel.empty()
    }
    all_consensus_reads = consensus_reads.mix(secondary_reads)

    // 5. consensus: spoa draft -> racon -> (medaka)
    SPOA_CONSENSUS(all_consensus_reads)
    MINIMAP2_ALIGN(all_consensus_reads.join(SPOA_CONSENSUS.out.draft))
    RACON(MINIMAP2_ALIGN.out.aligned)

    if (params.enable_medaka) {
        MEDAKA(all_consensus_reads.join(RACON.out.polished))
        consensus_ch = MEDAKA.out.consensus
    } else {
        consensus_ch = RACON.out.polished
    }

    // 6. QC: map ALL of the sample's usable reads (not just the dominant
    //    cluster) back to the consensus -- reads that don't map are the
    //    contamination/off-target signal, independent of the clustering. A second
    //    consensus is checked against its own cluster's reads instead: the rest
    //    of the sample is, by definition, not it
    MAPBACK(consensus_ch.join(usable_ch.mix(secondary_reads)))
    CONSENSUS_QC(MAPBACK.out.paf)

    // 7. one summary row per sample (including samples that got no consensus),
    //    a combined FASTA, and the QC html
    BUILD_REPORT(
        DOWNSAMPLE.out.counts.map { sample, tsv -> tsv }.collect(),
        purity_files,
        CONSENSUS_QC.out.qc.map { sample, tsv -> tsv }.collect().ifEmpty([]),
        consensus_ch.map { sample, fasta -> fasta }.collect().ifEmpty([])
    )

    emit:
    consensus = consensus_ch
    report    = BUILD_REPORT.out.summary
}
