#!/usr/bin/env nextflow

process FRAG_LONGER_100 {
    label 'process_medium'

    container = 'oras://community.wave.seqera.io/library/samtools:1.22--105e5e643c53f059'
    publishDir "bam.bt2", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.Longer100bp.bam"), path("*.Longer100bp.bam.bai"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -h ${bam} \\
    | awk 'substr(\$0,1,1)=="@" || (\$9 >= 100) || (\$9 <= -100)' \\
    | samtools view -@ ${task.cpus} -hb - > ${bam.baseName}.Longer100bp.bam

    samtools index -@ ${task.cpus} ${bam.baseName}.Longer100bp.bam
    """
}
