
library(readxl)
measurments <- read_excel("/rds/projects/c/croftap-labdata2/Chris/data_plot_RUNX1_paper/FigureS6A.xlsx")

summary_data <- measurments %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(percentage),
    SD = sd(percentage)
  )


y_max <- max(measurments$percentage, na.rm = TRUE)


summary_data$Condition <- factor(summary_data$Condition, levels = c("EV", "RUNX1a", "RUNX1c"))
measurments$Condition <- factor(measurments$Condition, levels = c("EV", "RUNX1a", "RUNX1c"))

summary_data %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill=c("red", "red", "red")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +        # Error bars
  geom_jitter(data = measurments, aes(x = Condition, y = percentage), width = 0.15,         # Points
              color = "black", size = 3.5, alpha = 0.8) +
   theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) 

