#!/bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 68:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear

module load bear-apps/2021b
module load SAMtools/1.15.1-GCC-11.2.0

cd /rds/projects/m/mahonyc-gata2wgs/Runx1_ChipSeq_GSE29180

samtools merge chipseqRUNX1_merged.bam unique_SRR443852_sorted.bam unique_SRR443853_sorted.bam unique_SRR443854_sorted.bam
samtools index chipseqRUNX1_merged.bam

