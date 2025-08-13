library(GenomicRanges)

RUNX1_peaks_homer <- read.delim("/rds/projects/m/mahonyc-runx1-bulk-seq-data/CHIPseq/RUNX1_peaks_homer.txt", header=FALSE)

ATAC_HIST_peaks_homer <- read.delim("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/Bulk_human_data/GSE128642_paper_ChipSeq/ATAC_HIST_peaks_homer.txt", header=FALSE)

gr1 <- GRanges(seqnames = RUNX1_peaks_homer$V1,
               ranges = IRanges(start = RUNX1_peaks_homer$V2, end = RUNX1_peaks_homer$V3))

gr2 <- GRanges(seqnames = ATAC_HIST_peaks_homer$V2,
               ranges = IRanges(start = ATAC_HIST_peaks_homer$V3, end = ATAC_HIST_peaks_homer$V4))


overlaps <- findOverlaps(gr1, gr2)

overlapping_gr1 <- gr1[queryHits(overlaps)]

overlapping_gr2 <- gr2[subjectHits(overlaps)]


overlapping_regions <- data.frame(
  chr1 = seqnames(overlapping_gr1),
  start1 = start(overlapping_gr1),
  end1 = end(overlapping_gr1),
  chr2 = seqnames(overlapping_gr2),
  start2 = start(overlapping_gr2),
  end2 = end(overlapping_gr2)
)



overlapping_regions[,c(1:3)] %>% write.table(file = "/rds/projects/m/mahonyc-runx1-bulk-seq-data/CHIPseq/overlap.csv", append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")


RUNX1_peaks_overlap_homer_anno <- read.delim("/rds/projects/m/mahonyc-runx1-bulk-seq-data/CHIPseq/RUNX1_peaks_overlap_homer_anno.txt")

unique(RUNX1_peaks_overlap_homer_anno$Gene.Name) %>% length()

RUNX1_peaks_overlap_homer_anno_f <- RUNX1_peaks_overlap_homer_anno %>% filter(grepl("promoter", Annotation))

unique(RUNX1_peaks_overlap_homer_anno_f$Gene.Name) %>% length()

library(enrichR)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")

enriched <- enrichr(unique(RUNX1_peaks_overlap_homer_anno$Gene.Name), dbs)

plotEnrich(enriched[[3]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")



enriched[[3]]$log_adj_p <- -log(enriched[[3]]$Adjusted.P.value)


enriched[[3]]$sig <- "NO"
enriched[[3]]$sig[enriched[[3]]$Adjusted.P.value > 0.05] <- "no"
enriched[[3]]$sig[enriched[[3]]$Adjusted.P.value < 0.05] <- "yes"

goterms_check <- enriched[[3]] %>% filter(log_adj_p > 1)


goterms_check <- enriched[[3]] %>% filter(Term %in% c("regulation of cartilage development (GO:0061035)", "regulation of cell-matrix adhesion (GO:0001952)", "positive regulation of cartilage development (GO:0061036)", "response to transforming growth factor beta (GO:0071559)", "regulation of interleukin-12 biosynthetic process (GO:0045075)", "regulation of fibroblast growth factor receptor signaling pathway (GO:0040036)", "positive regulation of cytokine production involved in immune response (GO:0002720)", "platelet-derived growth factor receptor signaling pathway (GO:0048008)", "substrate-dependent cell migration (GO:0006929)"))

plotEnrich(goterms_check, showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")

library(enrichplot)
dotplot(goterms_check)


enriched[[3]] %>%  ggplot(aes(x= Odds.Ratio, y= log_adj_p, label=Term)) + 
                  geom_point() +theme_ArchR()+
    geom_point(data = enriched[[3]] %>% filter(Term == "inflammatory response (GO:0006954)" ), color = "red")
colnames(enriched[[3]])


RUNX1_peaks_overlap_homer_anno$Annotation2 <- RUNX1_peaks_overlap_homer_anno$Annotation
library(splitstackshape)
RUNX1_peaks_overlap_homer_anno <- cSplit(RUNX1_peaks_overlap_homer_anno, splitCols = "Annotation2", sep=" ")


df1_pt <- table(RUNX1_peaks_overlap_homer_anno$Annotation2_1) %>% as.data.frame() 
df1_pt$pt <- df1_pt$Freq/sum(df1_pt$Freq)*100

df1_pt %>% ggplot( aes(x = "", y = pt, fill=Var1)) + geom_col(width=1, color=1) + coord_polar(theta = "y") +theme_ArchR()+ggtitle("RUNX1 ChIPseq peak annotations")

grn <- read_csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/grn.csv")

grn_RUNX1 <- grn %>% filter(Source=="RUNX1")

all_DEGs_human <-all_DEGs_human %>% filter(comparison == "EV_vs_R1C" & padj <0.05 & log2FoldChange < -0.01)


RUNX1_peaks_overlap_homer_anno_bulk <- RUNX1_peaks_overlap_homer_anno[RUNX1_peaks_overlap_homer_anno$Gene.Name %in% all_DEGs_human$gene_name,]


RUNX1_peaks_overlap_homer_anno_bulk$Gene.Name %>% unique() %>% length

library(enrichR)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")

enriched <- enrichr(unique(RUNX1_peaks_overlap_homer_anno_bulk$Gene.Name), dbs)

plotEnrich(enriched[[3]], showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")

goterms_check <- enriched[[3]] %>% filter(Term %in% c("positive regulation of cytokine production (GO:0001819)", "regulation of cytokine production (GO:0001817)", "lymphocyte activation (GO:0046649)",  "regulation of cell activation (GO:0050865)", "bone remodeling (GO:0046849)", "chemotaxis (GO:0006935)"))

plotEnrich(goterms_check, showTerms = 30, numChar = 40, y = "Count", orderBy = "P.value")
