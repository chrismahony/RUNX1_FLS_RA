setwd("/rds/projects/c/croftap-lab-data-2021/Paulynn/FORTESSA/MS RUNX1 KO STIA/21.05.25/SG/Sg Runx1 KO fib/SG Runx1 KOP ms fib")
library(readxl)
analysis_R <- read_excel("CM_analysis.xlsx")

analysis_R <- analysis_R %>% dplyr::filter(Percentage != 0.60)

summary_data <- analysis_R %>%
  dplyr::group_by(Condition, Cell) %>%
  dplyr::summarise(
    Mean = mean(Percentage),
    SD = sd(Percentage)
  )

for (i in 1:length(unique(summary_data$Cell))){
  
  analysis_R_f <- analysis_R %>% dplyr::filter(Cell == unique(summary_data$Cell)[i])
  
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
  

print(summary_data %>% dplyr::filter(Cell == unique(summary_data$Cell)[i]) %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.95, fill=c("lightgrey", "red")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +        # Error bars
  geom_jitter(data = analysis_R_f, aes(x = Condition, y = Percentage), width = 0.15,         # Points
              color = "black", size = 2.5, alpha = 0.8) +  # Points
    theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      )+
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(unique(analysis_R_f$Cell))
+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.5,  # space above highest point
           label = sig_label,
           size = 6))



print(summary_data %>% dplyr::filter(Cell == unique(summary_data$Cell)[i]) %>% 
  ggplot(aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.95, fill = "lightgrey")  +
  geom_jitter(data = analysis_R_f,
              aes(x = Condition, y = Percentage, fill = Condition),
              width = 0.15,
              shape = 21,            # Needed to allow fill + border
              color = "black",       # Border color
              size = 3.5,
              alpha = 1) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2)+
  theme_minimal() +
  theme_ArchR() +
  annotate("text", 
           x = 1.5,
           y = y_max + 0.001,
           label = sig_label,
           size = 6) + scale_fill_manual(values = c("TAM" = "#E41A1C", "Ctrl" = "#377EB8"))+
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(unique(analysis_R_f$Cell)))



}
