#!/usr/bin/env nextflow

process FRAG_SHORTER_100 {
    label 'process_medium'

    container = 'oras://community.wave.seqera.io/library/samtools:1.22--105e5e643c53f059'
    publishDir "results/bam.bt2", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.Shorter100bp.bam"), path("*.Shorter100bp.bam.bai"), emit: bam_bai

    script:
    """
    samtools view -@ ${task.cpus} -h ${bam} \\
    | awk 'substr(\$0,1,1)=="@" || (\$9 < 100 && \$9 > 0) || (\$9 > -100 && \$9 < 0)' \\
    | samtools view -@ ${task.cpus} -hb - > ${bam.baseName}.Shorter100bp.bam

    samtools index -@ ${task.cpus} ${bam.baseName}.Shorter100bp.bam
    """
}
