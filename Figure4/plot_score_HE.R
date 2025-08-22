library(readxl)
measurments <- read_excel("/rds/projects/c/croftap-labdata2/Chris/20250723_H&E/scoreing.xlsx")

summary_data <- measurments %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(`Total Synovitis`),
    SD = sd(`Total Synovitis`)
  )


y_max <- max(measurments$`Total Synovitis`, na.rm = TRUE)


t_res <- t.test(`Total Synovitis` ~ Condition, data = measurments)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "***",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)



summary_data$Condition <- factor(summary_data$Condition, levels = c("CTRL", "EXP"))
measurments$Condition <- factor(measurments$Condition, levels = c("CTRL", "EXP"))

summary_data %>% 
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("black", "red"), width = 0.6, fill = c("white")) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, color = c("black", "red")) +          # Error bars
  geom_jitter(data = measurments, aes(x = Condition, y = `Total Synovitis`, color=Condition), width = 0.15,       # Points
              size = 2.5, alpha = 1) +
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black")
  ) +
  scale_color_manual(values = c("CTRL" = "black", "EXP" = "red"))+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.1,  # space above highest point
           label = sig_label,
           size = 6)  +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+NoLegend()
