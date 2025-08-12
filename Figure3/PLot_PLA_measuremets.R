library(readxl)
cell_detections <- read_excel("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/RUNX1_inhibitor/PLA_25_6_2025/cell_detections.xlsx")

summary_data <- cell_detections %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(PLA_IF_cell),
    SD = sd(PLA_IF_cell)
  )

t_res <- t.test(PLA_IF_cell ~ Condition, data = cell_detections)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "***",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(cell_detections$PLA_IF_cell, na.rm = TRUE)


summary_data$Condition <- factor(summary_data$Condition, levels = c("DMSO", "CBFBI"))


summary_data %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("red", "black"), width = 0.55, fill=c("white"))  +        # Error bars
  geom_jitter(data = cell_detections, aes(x = Condition, y = PLA_IF_cell, color=Condition), width = 0.15, size = 2.5, alpha = 1) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = c("red", "black"))+
  scale_color_manual(values = c("DMSO" = "black", "CBFBI" = "red"))+
    # Remove gap between bar and x-axis
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black")
  )+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max+0.001,  # space above highest point
           label = sig_label,
           size = 6)+
  ggtitle("Cell/mm2")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+NoLegend()
