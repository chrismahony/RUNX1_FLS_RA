nonstim_merged_loops <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/nonstim_merged_loops.bedpe")

TNFa_merged_loops <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/TNFa_merged_loops.bedpe")

`8mix_merged_loops` <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/8mix_merged_loops.bedpe")

library("readr")
nonstim_merged_loops_runx1only<-nonstim_merged_loops[nonstim_merged_loops$chr1 == '21' & nonstim_merged_loops$x1 > 36160098 & nonstim_merged_loops$x1 < 36421595,]

write_delim(nonstim_merged_loops_runx1only, "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/nonstim_merged_loops_runx1only.bedpe", delim = "\t")

TNFa_merged_loops_runx1<-TNFa_merged_loops[TNFa_merged_loops$chr1 == '21' & TNFa_merged_loops$x1 > 36160098 & TNFa_merged_loops$x1 < 36421595,]

write_delim(TNFa_merged_loops_runx1, "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/TNFa_merged_loops_runx1.bedpe", delim = "\t")

`8mix_merged_loops_RUNX1`<-`8mix_merged_loops`[`8mix_merged_loops`$chr1 == '21' & `8mix_merged_loops`$x1 > 36160098 & `8mix_merged_loops`$x1 < 36421595,]

write_delim(`8mix_merged_loops_RUNX1`, "/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/hum0207.v1_Hi-C/hum0207.v1.HiC.v1/Hi-C/T8mix_merged_loops_RUNX1.bedpe", delim = "\t")
