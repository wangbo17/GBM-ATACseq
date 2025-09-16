#!/usr/bin/env nextflow

process BAM_TO_BIGWIG {
    label 'process_medium'
    container = 'oras://community.wave.seqera.io/library/deeptools:3.5.6--718d316362a5e588'

    publishDir "results/bigwig", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bigwig"), emit: bigwig

    script:
    """
    bamCoverage \\
        --bam ${bam} \\
        --outFileName ${bam.baseName}.bigwig \\
        --outFileFormat bigwig \\
        --scaleFactor 1 \\
        --binSize 1 \\
        --numberOfProcessors ${task.cpus} \\
        --normalizeUsing RPKM
    """
}
