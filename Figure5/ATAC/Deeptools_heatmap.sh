#!/bin/bash
#SBATCH -n 60
#SBATCH -N 1
#SBATCH --mem 180000
#SBATCH --time 48:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --account=croftap-stia-atac

#SRR5077625=H3K4me3 MEF
#SRR5077641=MEF_H3K27ac
#SRR5077728=Runx1 chipseq MEF

set -e

module purge; module load bluebear


module load deepTools/3.5.0-foss-2020a-Python-3.8.2

cd /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF

computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-R /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/ALL_PEAKS_728_peaks_heatmap.txt \
-S /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/ATACseq.bw \
--skipZeros \
-p 6 \
-o /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/ATACmatrix.gz \
--outFileSortedRegions /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/regionsATAC.bed \

plotHeatmap -m /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/ATACmatrix.gz \
      -out ATAC.png \

computeMatrix reference-point --referencePoint center \
-b 1000 -a 1000 \
-R /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/ALL_PEAKS_728_peaks_heatmap.txt \
-S /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/HACChipseq.bw \
--skipZeros \
-p 6 \
-o /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/HACChipseqmatrix.gz \
--outFileSortedRegions /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/regionsHACChipseq.bed \

plotHeatmap -m /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/GSE90893_Runx1_MEF/bigWig/HACChipseqmatrix.gz \
      -out HAC.png 
