library(readxl)
measurments <- read_excel("/rds/projects/c/croftap-labdata2/Chris/data_plot_RUNX1_paper/FigureS5E.xlsx")

summary_data <- measurments %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(`Live cells`),
    SD = sd(`Live cells`)
  )

t_res <- t.test(`Live cells` ~ Condition, data = measurments)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "***",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(measurments$`Live cells`, na.rm = TRUE)


summary_data$Condition <- factor(summary_data$Condition, levels = c("DMSO", "CBFβi"))
measurments$Condition <- factor(measurments$Condition, levels = c("DMSO", "CBFβi"))

summary_data %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("red", "black"), width = 0.55, fill="white") +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, color = c("red", "black")) +        # Error bars
  geom_jitter(data = measurments, aes(x = Condition, y = `Live cells`, color = Condition), width = 0.25,         # Points
               size = 3, alpha = 1)+
  scale_color_manual(values = c("DMSO" = "black", "CBFβi" = "red")) +
   theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black")
      ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+NoLegend()

