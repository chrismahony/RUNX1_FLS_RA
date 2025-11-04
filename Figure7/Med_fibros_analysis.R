
###HEatmap
library(readr)
grn <- read.csv("/rds/projects/m/mahonyc-cesar-data/fibro_analysis/scMEGA_repeat/df_grn2_new_FINAL.csv")

library(babelgene)
grn_Runx1 <- grn %>% filter(tf=="RUNX1")
grn_Runx1$gene <-   gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(grn_Runx1$gene), perl=TRUE)

human_genes <- orthologs(grn_Runx1$gene, species="mouse", human = FALSE)

med_fibros <- AddModuleScore(med_fibros, features=list(human_genes$human_symbol), name="RUNX1_grn")

Idents(med_fibros) <- 'Tissue'
med_fibros_no_lung <- subset(med_fibros, idents= levels(med_fibros)[-3])
med_fibros_no_lung <- med_fibros_no_lung %>% ScaleData()

Idents(med_fibros_no_lung)<-"Cluster_name"
Dotplot_heatmap<-DotPlot(med_fibros_no_lung, features = c("RUNX1", "RUNX1_grn1"))
Dotplot_heatmap_data <- Dotplot_heatmap[["data"]]

dotplot<-Dotplot_heatmap_data %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, border=T)



## Plot Corrolation Scatterplot with RUNX1 expression and inflam score

med_fibros_no_lung$SampleID %>% unique()
#extract avg. scaled expression for your gene in each cluster and sample
Idents(med_fibros_no_lung)<-"SampleID"
dotplot<-DotPlot(med_fibros_no_lung, features = "RUNX1_grn1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


#dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")
library(splitstackshape)
inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

dotplot_data$id_1 <- dotplot_data$id

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_1))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") #+
        #facet_wrap(~id_2)


ml = lm(Runx1~inflam, data = final_df)
summary(ml)$r.squared

library(ggpubr)

ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = -2, label.y = 1)



ggscatter(final_df, x = "Runx1", y = "inflam",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "black", fill = NA),
          shape = 21,             # Hollow circles
          color = "red",        # Border color
          fill = "white",         # Inside fill color (optional)
          size = 2 ,
           stroke = 1.6# Adjust size as needed
) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1)



### PLot RUNX1 GRN with inflam score 

Idents(med_fibros_no_lung)<-"SampleID"
dotplot<-DotPlot(med_fibros_no_lung, features = "RUNX1")
dotplot_data<-dotplot[["data"]]
dotplot_data <- subset(dotplot_data, select = c(avg.exp.scaled, id))
names(dotplot_data)[names(dotplot_data)=="avg.exp.scaled"] <- "Runx1"


#dotplot_data <- dotplot_data %>% cSplit(splitCols = "id", sep=".")

inflam_df <- paste(med_fibros$SampleID, med_fibros$InflamScore, sep=".") %>% unique() %>% as.data.frame() %>% cSplit(splitCols = ".", sep=".")

colnames(inflam_df) <- c("id_1", "inflam", "extra")

inflam_df <- inflam_df %>% replace(is.na(.), 0)

inflam_df$inflam <- paste(inflam_df$inflam, inflam_df$extra, sep=".")
inflam_df$inflam <- as.double(inflam_df$inflam)

dotplot_data$id_1 <- dotplot_data$id

final_df <- dotplot_data %>% 
  left_join(inflam_df, by="id_1")


#plot
#dotplot_data_synovium <- dotplot_data[dotplot_data$condition_2 == "Synovium",]
ggplot(final_df, aes(x = Runx1, y = inflam)) +
    geom_point(aes(color = factor(id_1))) +
    stat_smooth(method = "lm",
        col = "black",
        se = FALSE,
        size = 0.5)+theme_ArchR()+ theme (legend.position = "none") #+
        #facet_wrap(~id_2)


ml = lm(Runx1~inflam, data = final_df)
summary(ml)$r.squared

library(ggpubr)

ggscatter(final_df, x = "Runx1", y = "inflam",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   )+ stat_cor(method = "pearson", label.x = -1, label.y = 1.1)



ggscatter(final_df, x = "Runx1", y = "inflam",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "black", fill = NA),
          shape = 21,             # Hollow circles
          color = "red",        # Border color
          fill = "white",         # Inside fill color (optional)
          size = 2 ,
           stroke = 1.6# Adjust size as needed
) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
