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
  geom_bar(stat = "identity", color = c("black", "red"), width = 0.6, fill=c("white")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color=c("black", "red")) +        # Error bars
  geom_jitter(data = analysis_R_f, aes(x = Condition, y = Percentage, color=Condition), width = 0.15, size = 2.5, alpha = 1) +  # Points
    theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black")
      )+
  scale_color_manual(values = c("Ctrl" = "black", "TAM" = "red")) +
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(unique(analysis_R_f$Cell))
+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.5,  # space above highest point
           label = sig_label,
           size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+NoLegend())
