library(readxl)
measurments <- read_excel("/rds/projects/c/croftap-labdata2/Chris/data_plot_RUNX1_paper/Figure4H.xlsx")

summary_data <- measurments %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(`Final score`),
    SD = sd(`Final score`)
  )


y_max <- max(measurments$`Final score`, na.rm = TRUE)


t_res <- wilcox.test(`Final score` ~ Condition, data = measurments)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "***",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)



summary_data$Condition <- factor(summary_data$Condition, levels = c("CTRL", "RUNX1_KO"))
measurments$Condition <- factor(measurments$Condition, levels = c("CTRL", "RUNX1_KO"))

summary_data %>% 
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("black", "red"), width = 0.6, fill = c("white")) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, size = 0.8, color = c("black", "red")) +          # Error bars
  geom_jitter(data = measurments, aes(x = Condition, y = `Final score`, color=Condition), width = 0.15,       # Points
              size = 3, alpha = 1) +
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = c("CTRL" = "black", "RUNX1_KO" = "red"))+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.1,  # space above highest point
           label = sig_label,
           size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+NoLegend()
