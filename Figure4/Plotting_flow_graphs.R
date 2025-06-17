setwd("/rds/projects/c/croftap-lab-data-2021/Paulynn/FORTESSA/MS RUNX1 KO STIA/08.04.25/Ms Runx1 flox STIA")
library(readxl)
analysis_R <- read_excel("analysis_R.xlsx")


summary_data <- analysis_R %>%
  dplyr::group_by(Condition, Cell) %>%
  dplyr::summarise(
    Mean = mean(Percentage),
    SD = sd(Percentage)
  )

for (i in 1:length(unique(summary_data$Cell))){
  
  analysis_R_f <- analysis_R %>% filter(Cell == unique(summary_data$Cell)[i])
  
      # Perform t-test
t_res <- t.test(Percentage ~ Condition, data = analysis_R_f)

# Get p-value and significance label
p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(analysis_R_f$Percentage, na.rm = TRUE)
  

print(summary_data %>% filter(Cell == unique(summary_data$Cell)[i]) %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.95, fill=c("lightgrey", "red")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +        # Error bars
  geom_jitter(data = analysis_R_f, aes(x = Condition, y = Percentage), width = 0.15,         # Points
              color = "black", size = 2.5, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(unique(analysis_R_f$Cell))
+theme_ArchR()+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.5,  # space above highest point
           label = sig_label,
           size = 6))

}
