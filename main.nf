#!/usr/bin/env nextflow
/*
 * Single-amplicon ONT consensus pipeline.
 * One barcode gene per sample (one sample per barcode) -> one polished
 * consensus sequence per sample, with QC.
 * Sibling of edna-ont-nf / 16S-nf (mixed amplicons) and ulana-nf (genomes);
 * reminiscent of epi2me-labs/wf-amplicon.
 */

nextflow.enable.dsl = 2

include { SINGLE_AMPLICON } from './workflows/single_amplicon.nf'
include { MERGE_FASTQ     } from './modules/merge_fastq.nf'
include { CHOPPER         } from './modules/chopper.nf'
include { READ_STATS; READ_STATS_REPORT } from './modules/read_stats.nf'

// every param default lives in nextflow.config (see the note there)

def helpMessage() {
    log.info """
    Single-amplicon ONT consensus pipeline
    ---------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results -profile docker

    Required:
      --input       CSV: sample,fastq (fastq may be a glob over a barcode's part-files)

    Key optional:
      --fwd_primer / --rev_primer   primer sequences (5'->3' as ordered); both needed to trim
      --min_len / --max_len / --min_qual   chopper filtering (default ${params.min_len} / ${params.max_len ?: 'none'} / ${params.min_qual}) -- set for your amplicon
      --downsample_reads   reads used for the consensus, 0 = all (default ${params.downsample_reads})
      --min_reads          samples with fewer reads are reported as FAIL, not assembled (default ${params.min_reads})
      --enable_purity_check  build the consensus from the dominant read cluster only (default ${params.enable_purity_check})
      --cluster_id         vsearch identity for the purity check (default ${params.cluster_id})
      --enable_secondary_consensus  also build a consensus from the second-largest read cluster (default ${params.enable_secondary_consensus}), if it holds at least --min_secondary_frac of the reads (default ${params.min_secondary_frac})
      --enable_medaka      medaka polish (default ${params.enable_medaka}); --medaka_model must match your basecaller (default ${params.medaka_model})
      --min_depth / --min_primary_frac / --min_dominant_frac   QC thresholds (default ${params.min_depth} / ${params.min_primary_frac} / ${params.min_dominant_frac})
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob matching zero files, so check explicitly
            if (row.sample.endsWith(params.secondary_suffix)) {
                error "Samplesheet row '${row.sample}': sample names can't end in '${params.secondary_suffix}', which the pipeline uses for a sample's second consensus"
            }
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

workflow {
    if (params.help || !params.input) {
        helpMessage()
        exit 0
    }

    SINGLE_AMPLICON(samplesheetToChannel(params.input))
}

// Read QC only: merge -> chopper -> before/after length & Q-score summary, no
// consensus. Handy for picking --min_len/--max_len/--min_qual before a full run:
//   nextflow run main.nf -entry READ_QC_ONLY --input samplesheet.csv -profile docker
workflow READ_QC_ONLY {
    if (!params.input) {
        log.error "READ_QC_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
    CHOPPER(MERGE_FASTQ.out.merged)
    READ_STATS(MERGE_FASTQ.out.merged.join(CHOPPER.out.filtered))
    READ_STATS_REPORT(
        READ_STATS.out.stats.map { sample, stats, hist -> stats }.collect(),
        READ_STATS.out.stats.map { sample, stats, hist -> hist }.collect()
    )
}
