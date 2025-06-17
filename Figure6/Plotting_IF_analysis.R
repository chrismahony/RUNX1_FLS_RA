measurments <- read_excel("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/RUNX1_IF/images_quantifiedmeasurments.xlsx")

measurments$cell_image <- paste(measurments$`Cell number`, measurments$Image, sep="_")

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


summary_data$Condition <- factor(summary_data$Condition, levels = c("Unstimulated", "DLL4 treated"))
measurments$Condition <- factor(measurments$Condition, levels = c("Unstimulated", "DLL4 treated"))

summary_data %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.95, fill=c("red", "lightgrey")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +        # Error bars
  geom_jitter(data = measurments, aes(x = Condition, y = `Mean intensity`), width = 0.15,         # Points
              color = "black", size = 2.5, alpha = 0.8) +
  theme_minimal() +
theme_ArchR()+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.5,  # space above highest point
           label = sig_label,
           size = 6)

