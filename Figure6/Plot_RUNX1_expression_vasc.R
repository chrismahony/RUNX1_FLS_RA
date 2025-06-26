library(readxl)
fnames <- list.files(path = "/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/Organoids/Synovium images", pattern = "qpath_image")

qpath_values <- list()

for (i in 1:length(fnames)){
  qpath_values[[i]] <- read_excel(fnames[[i]])
qpath_values[[i]] <- qpath_values[[i]] %>% dplyr::select(c(Image, Classification, `Distance to annotation with Vascular µm`, `Cell: ChS2-T2 mean`))
}

names(qpath_values) <- gsub(pattern=".xlsx",x=fnames, replacement="")

library(data.table)

# Synovium test 4.czi not used for distance analysies, many bg cells which were hard to remove.

rbindlist(qpath_values) %>% 
  filter(Classification %in% c("SL")) %>%
  filter(Image %in% c("Synovium test 1.czi", "Synovium test 2.czi", "Synovium test 3.czi")) %>% 
ggplot(aes(x=`Cell: ChS2-T2 mean`, y=`Distance to annotation with Vascular µm`))+
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  )+
  geom_smooth(color = "red")+xlim(0,1500)#+facet_wrap("Image")
