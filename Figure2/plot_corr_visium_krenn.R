setwd("/rds/projects/c/croftap-runx1data01/analysis_of_public_data_sets/visium_RUNX1_corr")
library(readxl)
data <- read_excel("data.xlsx")

ggscatter(data, x = "RUNX1 expression", y = "Krenn_Global",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "black", fill = NA),
          shape = 21,             # Hollow circles
          color = "red",        # Border color
          fill = "white",         # Inside fill color (optional)
          size = 2.5 ,
           stroke = 1.6# Adjust size as needed
) +
  stat_cor(method = "pearson", label.x = 1, label.y = 1)
