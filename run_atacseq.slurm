#!/bin/bash
#SBATCH --job-name=atacseq
#SBATCH --output=atacseq.out
#SBATCH --error=atacseq.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G

conda activate nf_env
nextflow run atacseq.nf -resume