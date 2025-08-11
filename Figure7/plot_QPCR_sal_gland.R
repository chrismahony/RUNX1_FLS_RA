library(readxl)
qPCR_data <- read_excel("/rds/projects/c/croftap-labdata2/Chris/data_plot_RUNX1_paper/qPCR_data.xlsx")

qPCR_data$Condition <- c(rep("ctrl",2), rep("-Runx1",3))
                         

summary_data <- qPCR_data %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(`2_neg_d_CTx100000`),
    SD = sd(`2_neg_d_CTx100000`)
  )


y_max <- max(qPCR_data$`2_neg_d_CTx100000`, na.rm = TRUE)


t_res <- wilcox.test(`2_neg_d_CTx100000` ~ Condition, data = qPCR_data)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "***",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

summary_data$Condition <- factor(summary_data$Condition, levels = c("ctrl", "-Runx1"))
qPCR_data$Condition <- factor(qPCR_data$Condition, levels = c("ctrl", "-Runx1"))


summary_data %>% 
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill = c("red", "grey")) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +          # Error bars
  geom_jitter(data = qPCR_data, aes(x = Condition, y = `2_neg_d_CTx100000`), width = 0.15,       # Points
              color = "black", size = 3.5, alpha = 0.8) +
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black"),# keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  ) +
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.1,  # space above highest point
           label = sig_label,
           size = 6) 
