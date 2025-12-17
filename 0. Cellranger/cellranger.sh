#!/bin/bash
#SBATCH --job-name=wwy_cnvkit
#SBATCH --error=%x.%j.err
#SBATCH --output=%x.%j.out
#SBATCH --partition=cpu1,cpu2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64000
outdir = /data/person/g5/wangwy/
cd outdir 
cellranger count \
--id run_pbmc_1k_v3_fastqs \
--create-bam true \
--transcriptome refdata-gex-GRCh38-2024-A/ \
--fastqs pbmc_1k_v3_fastqs/ \
--nosecondary \
--localcores 16 \
--localmem 64
