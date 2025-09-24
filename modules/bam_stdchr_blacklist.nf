#!/usr/bin/env nextflow

process BAM_STDCHR_BLACKLIST {
    label 'process_low'

    container = 'oras://community.wave.seqera.io/library/bedtools_samtools:c3dc2897ee2af45d'
    publishDir "results/bedtools", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fasta
    path blacklist

    output:
    tuple val(meta), path("*.filtered.dedup.clean.bam"), emit: bam

    script:
    def sample_id = meta.sample_id

    """
    samtools faidx ${fasta}
    cut -f1,2 ${fasta}.fai > chrom.sizes

    awk 'BEGIN{OFS="\\t"} \$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)\$/ {print \$1, 0, \$2}' chrom.sizes > stdchr.bed

    sort -k1,1 -k2,2n stdchr.bed > stdchr.sorted.bed
    sort -k1,1 -k2,2n ${blacklist} > blacklist.sorted.bed || true
    
    bedtools subtract -a stdchr.sorted.bed -b blacklist.sorted.bed \
      | bedtools merge -i - > allowed.bed

    bedtools intersect -abam ${bam} -b allowed.bed > ${sample_id}.filtered.dedup.clean.bam
    """
}
