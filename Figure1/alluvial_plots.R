
library(ggalluvial)
library(ggplot2)
library(Seurat)
library(ggforce)



Idents(stia2021_rna) <- 'cluster.name'
cluster.name <- Idents(stia2021_rna)  # Cluster identities
condition <- stia2021_rna$condition  # Condition information (replace 'condition' with your metadata column name)


data <- stia2021_rna@meta.data %>%
  group_by(cluster.name,condition) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )



data_labels <- tibble(
    group = c(
      rep('cluster.name', length(cluster.name)),
      rep('condition', length(condition))
    )
 ) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )


cluster_colors <- ArchR::paletteDiscrete(unique(data$cluster.name))

condition_colors <- c(
  "control" = "white",      # Bold orange-red, contrasting with fibroblast colors
  "initiation" = "white",   # Rich purple, stands out against blues
  "peak" = "white",         # Strong blue, a nice contrast to red and green fibroblasts
  "persistent" = "white",   # Bright green, stands out against yellow and orange
  "resolved" = "white",     # Yellow-orange, will differentiate well from purple
  "resolving" = "white"     # Strong red, contrasts well with blue and green fibroblasts
)

condition_colors <- c(
  "control" = "#FF5733",      # Bold orange-red, contrasting with fibroblast colors
  "initiation" = "#9B59B6",   # Rich purple, stands out against blues
  "peak" = "#1F77B4",         # Strong blue, a nice contrast to red and green fibroblasts
  "persistent" = "#2CA02C",   # Bright green, stands out against yellow and orange
  "resolved" = "#F39C12",     # Yellow-orange, will differentiate well from purple
  "resolving" = "#E74C3C"     # Strong red, contrasts well with blue and green fibroblasts
)

all_col <- c(cluster_colors, condition_colors)

ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = cluster.name), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +theme_ArchR()+
    scale_fill_manual(name = "Cluster", values = all_col)
            fontface = "bold") ) # Optional: make the text bold
  
  



Idents(stia2021_rna) <- 'cluster.name'
cluster.name <- Idents(stia2021_rna)  # Cluster identities
condition <- stia2021_rna$experiment  # Condition information (replace 'condition' with your metadata column name)


data <- stia2021_rna@meta.data %>%
  group_by(cluster.name,experiment) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )



data_labels <- tibble(
    group = c(
      rep('cluster.name', length(cluster.name)),
      rep('condition', length(condition))
    )
 ) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )


cluster_colors <- ArchR::paletteDiscrete(unique(data$cluster.name))

condition_colors <- c(
  "stia2017" = "white",      
  "stia2021" = "white",   
  "stia2018" = "white"        
  
)



all_col <- c(cluster_colors, condition_colors)

ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = cluster.name), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +theme_ArchR()+
    scale_fill_manual(name = "Cluster", values = all_col)
   


condition_colors <- c(
  "stia2017" = "#FF5733",      # Bold orange-red, contrasting with fibroblast colors
  "stia2021" = "#9B59B6",   # Rich purple, stands out against blues
  "stia2018" = "#1F77B4"         # Strong blue, a nice contrast to red and green fibroblasts
  
)

all_col <- c(cluster_colors, condition_colors)


ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = cluster.name), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +theme_ArchR()+
    scale_fill_manual(name = "Cluster", values = all_col)
  



Idents(stia2021_rna) <- 'cluster.name'
cluster.name <- Idents(stia2021_rna)  # Cluster identities
condition <- stia2021_rna$experiment  # Condition information (replace 'condition' with your metadata column name)


data <- stia2021_rna@meta.data %>%
  group_by(condition,experiment) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(x, levels = unique(x)),
    y = factor(y, levels = unique(y))
  )



data_labels <- tibble(
    group = c(
      rep('cluster.name', length(cluster.name)),
      rep('condition', length(condition))
    )
 ) %>%
  mutate(
    hjust = ifelse(group == 'sample', 1, 0),
    nudge_x = ifelse(group == 'sample', -0.1, 0.1)
  )


ggplot(data, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = condition), alpha = 0.75, axis.width = 0.15) +
  geom_parallel_sets_axes(aes(fill = y), color = 'black', axis.width = 0.1) +theme_ArchR()

