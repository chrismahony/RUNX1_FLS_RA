#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 68:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear

module load deepTools/3.5.0-foss-2020a-Python-3.8.2

cd /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF

#more details https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html

bamCoverage -b ./nomt_unique_SRR5077744_sorted.bam  \
-o ./bigWig/ATACseq.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> ./logs/Nanog_rep2_bamCoverage.log \

bamCoverage -b ./unique_SRR5077728_sorted.bam \
-o ./bigWig/Runx1Chipseq.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> ./logs/Nanog_rep2_bamCoverage.log  \

bamCoverage -b ./unique_SRR5077641_sorted.bam \
-o ./bigWig/HACChipseq.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 6 2> ./logs/Nanog_rep2_bamCoverage.log
