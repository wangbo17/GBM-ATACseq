#!/usr/bin/env nextflow

// Module INCLUDE statements
include { FASTQC } from './modules/fastqc.nf'
include { TRIM_GALORE } from './modules/trim_galore.nf'
include { BOWTIE2_BUILD } from './modules/bowtie2_build.nf'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align.nf'

include { BAM_PROPER_FILTER } from './modules/bam_proper_filter.nf'
include { PICARD_MARKDUPLICATES } from './modules/picard_markduplicates.nf'
include { BAM_MAPQ20_FILTER } from './modules/bam_mapq20_filter.nf'

include { SAMTOOLS_INDEX } from './modules/samtools_index.nf'
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf'
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats.nf'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat.nf'

include { FRAG_SHORTER_100 } from './modules/frag_shorter_100.nf'
include { FRAG_LONGER_100 } from './modules/frag_longer_100.nf'

include { BAM_TO_BIGWIG as BAM_TO_BIGWIG_SHORT } from './modules/bam_to_bigwig.nf'
include { BAM_TO_BIGWIG as BAM_TO_BIGWIG_LONG  } from './modules/bam_to_bigwig.nf'

include { MULTIQC } from './modules/multiqc.nf'

/*
 * Pipeline parameters
 */

// Primary input
params.input_csv = "data/samplesheet_test.csv"
params.fasta = "test-datasets/reference/genome.fa"
params.gtf = "test-datasets/reference/genes.gtf"
params.report_id = "atacseq"

log.info """\
==================================================================================================
\033[36m       ____ ____  __  __      _  _____  _    ____     ____             
      / ___| __ )|  \\/  |    / \\|_   _|/ \\  / ___|   / ___|  ___  __ _ 
     | |  _|  _ \\| |\\/| |   / _ \\ | | / _ \\| |   ____\\___ \\ / _ \\/ _` |
     | |_| | |_) | |  | |  / ___ \\| |/ ___ \\ |__|_____|__) |  __/ (_| |
      \\____|____/|_|  |_| /_/   \\_\\_/_/   \\_\\____|   |____/ \\___|\\__, |
                                                                    |_|\033[0m
                                                                                     
==================================================================================================

  \033[1;37m▸ Key Parameters:\033[0m

    ▹ \033[33mExecuted By\033[0m       : ${System.getProperty('user.name')}
    ▹ \033[33mRun Date\033[0m          : ${workflow.start.format("yyyy-MM-dd HH:mm 'UTC'")}

                                                                Author: Bo Wang | Version: Beta
==================================================================================================
"""
.stripIndent(true)

workflow {

    // INPUT CHANNEL CREATION

    read_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header:true)
        .map { row ->
            def meta = [
                sample_id: row.sample_id
            ]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }

    // PREPROCESSING AND ALIGNMENT

    FASTQC(read_ch)

    TRIM_GALORE(read_ch)

    BOWTIE2_BUILD(file(params.fasta))

    BOWTIE2_ALIGN(TRIM_GALORE.out.trimmed_reads, BOWTIE2_BUILD.out.bt2_index)

    // DUPLICATE REMOVAL AND FILTERING

    BAM_PROPER_FILTER(BOWTIE2_ALIGN.out.bam)

    PICARD_MARKDUPLICATES(BAM_PROPER_FILTER.out.bam, file(params.fasta))
    
    BAM_MAPQ20_FILTER(PICARD_MARKDUPLICATES.out.bam)

    // BAM INDEXING & QC STATISTICS

    SAMTOOLS_INDEX(BAM_MAPQ20_FILTER.out.bam)

    SAMTOOLS_STATS(SAMTOOLS_INDEX.out.bam_bai, file(params.fasta))

    SAMTOOLS_IDXSTATS(SAMTOOLS_INDEX.out.bam_bai)

    SAMTOOLS_FLAGSTAT(SAMTOOLS_INDEX.out.bam_bai)

    // FRAGMENT LENGTH SELECTION

    FRAG_SHORTER_100(SAMTOOLS_INDEX.out.bam_bai)

    FRAG_LONGER_100(SAMTOOLS_INDEX.out.bam_bai)

    BAM_TO_BIGWIG_SHORT(FRAG_SHORTER_100.out.bam)

    BAM_TO_BIGWIG_LONG(FRAG_LONGER_100.out.bam)

    // SUMMARY REPORT GENERATION

    MULTIQC(
        FASTQC.out.zip.mix(
            FASTQC.out.html,
            TRIM_GALORE.out.trimming_reports,
            TRIM_GALORE.out.fastqc_reports_1,
            TRIM_GALORE.out.fastqc_reports_2,
            BOWTIE2_ALIGN.out.log,
            PICARD_MARKDUPLICATES.out.metrics,
            SAMTOOLS_STATS.out.stats,
            SAMTOOLS_IDXSTATS.out.idxstats,
            SAMTOOLS_FLAGSTAT.out.flagstat
        ).collect(),
        params.report_id
    )
}