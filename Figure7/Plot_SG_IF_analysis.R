setwd("/rds/projects/c/croftap-labdata2/Chris/SG_slidescanner/stellaris/analysis")
library(readxl)
analysis <- read_excel("analysis.xlsx")

analysis$aggr_mm <- analysis$`aggr/total area` * 1000


summary_data <- analysis %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(aggr_mm),
    SD = sd(aggr_mm))
  


t_res <- t.test(aggr_mm ~ Condition, data = analysis)

# Get p-value and significance label
p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "****",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(analysis$aggr_mm, na.rm = TRUE)
  

summary_data %>%
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.95, 
           fill = c("lightgrey", "red")) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +  # Error bars
  geom_jitter(data = analysis, 
              aes(x = Condition, y = aggr_mm), 
              width = 0.15, color = "black", size = 2.5, alpha = 0.8) +  # Points
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  ) +
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value") +
  annotate("text", 
           x = 1.5, 
           y = max(summary_data$Mean + summary_data$SD) + 0,  # dynamically compute max y
           label = sig_label,
           size = 6)

