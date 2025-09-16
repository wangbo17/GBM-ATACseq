#!/usr/bin/env nextflow

process MACS2_CALLPEAK {
    label 'process_medium'
    container = 'oras://community.wave.seqera.io/library/macs2:2.2.7.1--189209832f83d10b'
    publishDir "results/macs2", mode: 'copy'

    input:
    tuple val(meta), path(exp_bam), path(exp_bai)

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peak
    tuple val(meta), path("*.xls"), emit: xls

    script:
    def sample_id = meta.sample_id

    """
    macs2 callpeak \\
        --treatment ${exp_bam} \\
        --format BAMPE \\
        --gsize hs \\
        --keep-dup all \\
        --bdg \\
        --SPMR \\
        --nomodel \\
        --shift -75 \\
        --extsize 150 \\
        --name ${sample_id}_ATAC
    """
}
