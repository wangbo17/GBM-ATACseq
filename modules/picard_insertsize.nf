#!/usr/bin/env nextflow

process PICARD_INSERTSIZE {
    label 'process_medium'

    container = 'oras://community.wave.seqera.io/library/picard_samtools:33d1f34d6faf154e'
    publishDir "results/picard/markduplicates", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.InsertSize.txt"), emit: metrics
    tuple val(meta), path("*.InsertSize.pdf"), emit: plot

    script:
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072

    def sample_id = meta.sample_id
    
    """
    picard \\
        -Xmx${avail_mem}M \\
        CollectInsertSizeMetrics \\
        --INPUT ${bam} \\
        --OUTPUT ${sample_id}.InsertSize.txt \\
        --HISTOGRAM_FILE ${sample_id}.InsertSize.pdf \\
        --METRIC_ACCUMULATION_LEVEL ALL_READS
    """
}