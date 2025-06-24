setwd("/rds/projects/c/croftap-runx1data01/luciferase_analysis")
library(readxl)
data <- read_excel("data.xlsx")

data$Condition_DOX <- paste(data$Condition, data$DOX, sep="_")

summary_data <- data %>%
  dplyr::group_by(Condition_DOX, Gene) %>%
  dplyr::summarise(
    Mean = mean(`Normalised Luciferase Activity`),
    SD = sd(`Normalised Luciferase Activity`)
  )


sig_df <- list()

for (i in 1:length(unique(summary_data$Gene))){
  
  analysis_R_f <- data %>% filter(Gene == unique(data$Gene)[i])
  
  
  aov_res <- aov(`Normalised Luciferase Activity` ~ Condition_DOX, data = analysis_R_f)
tukey_res <- TukeyHSD(aov_res)

tukey_df <- as.data.frame(tukey_res$Condition)
tukey_df$Comparison <- rownames(tukey_df)

# Split into group1 and group2 for plotting
tukey_df <- tukey_df %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(p.adj.signif = case_when(
    `p adj` < 0.001 ~ "****",
    `p adj` < 0.001 ~ "***",
    `p adj` < 0.01  ~ "**",
    `p adj` < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))


sig_df[[i]] <- tukey_df

y_max <- max(analysis_R_f$`Normalised Luciferase Activity`, na.rm = TRUE)+0.5
  

print(
  summary_data %>%
  filter(Gene == unique(summary_data$Gene)[i]) %>%
  ggplot(aes(x = Condition_DOX, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, 
           fill = c("lightgrey", "lightgrey", "red", "red")) +  # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +  # Error bars
  geom_jitter(data = analysis_R_f, 
              aes(x = Condition_DOX, y = `Normalised Luciferase Activity`), 
              width = 0.15, color = "black", size = 3.5, alpha = 0.8) +  # Remove gap between bar and x-axis
  theme(
    panel.background = element_blank(),       # remove inner panel background
    plot.background = element_blank(),        # remove outer background
    panel.grid.major = element_blank(),       # remove major grid lines
    panel.grid.minor = element_blank(),       # remove minor grid lines
    axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
  ) +
    ggtitle(unique(analysis_R_f$Gene)) 
)

}


names(sig_df) <- unique(summary_data$Gene)
