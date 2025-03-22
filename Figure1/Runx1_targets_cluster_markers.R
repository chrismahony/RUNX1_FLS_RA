df.grn2_RUNX1 <- df.grn2 %>% filter(tf == "RUNX1")

df.grn2_RUNX1$gene <-   gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(df.grn2_RUNX1$gene), perl=TRUE)

all_markers_f <- all_markers %>% filter(p_val_adj < 0.05)

all_markers_f_RUNX1_grn <- all_markers_f[all_markers_f$gene %in% df.grn2_RUNX1$gene,]


df_table <- data.frame(RUNX1_targets=table(all_markers_f_RUNX1_grn$cluster), all_genes=table(all_markers_f$cluster))

df_table$pct <- (df_table$RUNX1_targets.Freq/df_table$all_genes.Freq*100)

colours <- list('cluster' = ArchR::paletteDiscrete(stia2021_rna@meta.data[, "cluster.name"]),'condition'= ArchR::paletteDiscrete(stia2021_rna@meta.data[, "condition"]))


cols <- ArchR::paletteDiscrete(stia2021_rna@meta.data[, "cluster.name"]) %>% as.data.frame()

df_table_f <- df_table[df_table$RUNX1_targets.Var1 != "fibroblast__Runx2_Bglap",]
df_table_f <- df_table_f[df_table_f$RUNX1_targets.Var1 != "fibroblast_mural_Rgs5",]
df_table_f <- df_table_f[df_table_f$RUNX1_targets.Var1 != "fibroblast__Clu",]

p1 <- ggplot(df_table_f, aes(x = pct, y = reorder(RUNX1_targets.Var1, pct), fill = RUNX1_targets.Var1)) +
  geom_col() +
  scale_fill_manual(
                    values = cols$.) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "off")+theme_ArchR()+NoLegend()

unique(p1[["data"]][["RUNX1_targets.Var1"]])

levels(stia2021_rna) <- c("fibroblast_lining_F13a1_Col22a1"   ,  "fibroblast_sublining_Pi16"   ,"fibroblast__Crabp1_Col23a1","fibroblast_sublining_Serpina3c_C3"  ,  "fibroblast_sublining_Ccl11"        ,"fibroblast_sublining_Fmo2"   , "fibroblast__Chodl" ,"fibroblast_sublining_lining_Ccl7_Ccl2", "fibroblast_sublining_Sfrp1_Cfb"     ,   
 "fibroblast_sublining_C1qtnf3_Col8a1",                               
  "fibroblast__Runx2_Bglap", "fibroblast_mural_Rgs5", "fibroblast__Clu")




p2 <- DotPlot(stia2021_rna, features="Runx1", idents = unique(p1[["data"]][["RUNX1_targets.Var1"]]))+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10), text = element_text(size=10))

plot_grid(p1, p2)



ggplot(df_table_f, aes(x=RUNX1_targets.Freq, y=pct)) + 
  geom_point(size=2)+theme_ArchR()+geom_vline(xintercept = mean(df_table_f$RUNX1_targets.Freq), linetype = "dashed", color = "red") +  # Vertical red dashed line
  geom_hline(yintercept = mean(df_table_f$pct), linetype = "dashed", color = "red")+
  geom_text_repel(aes(label = df_table_f$RUNX1_targets.Var1
),  size = 3, color = "black") 


ggplot(df_table_f, aes(x=RUNX1_targets.Freq, y=pct)) + 
  geom_point(size=2)+theme_ArchR()+geom_vline(xintercept = mean(df_table_f$RUNX1_targets.Freq), linetype = "dashed", color = "red") +  # Vertical red dashed line
  geom_hline(yintercept = mean(df_table_f$pct), linetype = "dashed", color = "red")
