library(readxl)
measurments <- read_excel("//its-rds.bham.ac.uk/rdsprojects/c/croftap-stia-atac/CM_multiome/Functional_validation/RUNX1_IF/images_quantified/measurments_TNFa.xlsx")



measurments$cell_image <- paste(measurments$`Cell number`, measurments$Image, sep="_")

library(dplyr)
library(ggplot2)
summary_data <- measurments %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(`Mean intensity`),
    SD = sd(`Mean intensity`)
  )

t_res <- t.test(`Mean intensity` ~ Condition, data = measurments)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(measurments$`Mean intensity`, na.rm = TRUE)


summary_data$Condition <- factor(summary_data$Condition, levels = c("Unstimulated", "TNFa"))
measurments$Condition <- factor(measurments$Condition, levels = c("Unstimulated", "TNFa"))

summary_data %>%
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("red", "black"), fill = "white", width = 0.65) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, color = c("red", "black")) +  # Error bars
  geom_jitter(
    data = measurments, 
    aes(x = Condition, y = `Mean intensity`, color = Condition), 
    width = 0.25, size = 3, alpha = 1
  ) +
  scale_color_manual(values = c("Unstimulated" = "black", "TNFa" = "red")) +
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
  annotate(
    "text", 
    x = 1.5, 
    y = max(summary_data$Mean + summary_data$SD), 
    label = sig_label,
    size = 6
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))


summary_data %>%
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = c("red", "black"), fill = "white", width = 0.65) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.3, color = c("red", "black")) +  # Error bars
  geom_jitter(
    data = measurments, 
    aes(x = Condition, y = `Mean intensity`, color = Condition), 
    width = 0.25, size = 3, alpha = 1
  ) +
  scale_color_manual(values = c("Unstimulated" = "black", "TNFa" = "red")) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
