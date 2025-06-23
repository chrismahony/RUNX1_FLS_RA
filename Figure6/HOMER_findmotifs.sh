#!/bin/bash
#SBATCH --time 99:0:0
#SBATCH --qos castles
#SBATCH --mail-type ALL
#SBATCH --mem 300G
#SBATCH --ntasks 40
#SBATCH --account=croftap-stia-atac


set -e

module purge; module load bluebear

cd /rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1.ChIP.v1

PATH=/rds/projects/m/mahonyc-gata2-mutated-mds-aml/scATACseqAugust2020_analysis/homer/bin:${PATH}

findMotifsGenome.pl RA_all_homer.txt hg19 \
ALL_vs_NS -bg RA_ns_homer.txt  -bits -size given

findMotifsGenome.pl RA_TNFa_homer.txt hg19 \
TNFa_vs_NS -bg RA_ns_homer.txt  -bits -size given

