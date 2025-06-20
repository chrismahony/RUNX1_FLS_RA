setwd("/rds/projects/c/croftap-lab-data-2021/Paulynn/FORTESSA/RUNX1 OE Co-culture/01.05..25")
library(readxl)
scores <- read_excel("data_R.xlsx")


summary_data <- scores %>%
  dplyr::group_by(Group, Condition, Cell) %>%
  dplyr::summarise(
    Mean = mean(FC),
    SD = sd(FC)
  )


  summary_data$group_condition <- paste(summary_data$Group, summary_data$Cell)

scores$group_condition <- paste(scores$Group, scores$Cell)


for (i in 1:length(unique(summary_data$group_condition))){
  
  analysis_R_f <- scores %>% filter(group_condition == unique(summary_data$group_condition)[i])
  analysis_R_f$group_condition <- paste(analysis_R_f$Group, analysis_R_f$Cell)

      # Perform t-test
u_res <- wilcox.test(FC ~ Condition, data = analysis_R_f)

# Get p-value and significance label
p_val <- u_res$p.value
sig_label <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(analysis_R_f$FC, na.rm = TRUE)
  

print(summary_data %>% filter(group_condition == unique(summary_data$group_condition)[i]) %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill=c("lightgrey", "red")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +        # Error bars
  geom_jitter(data = analysis_R_f, aes(x = Condition, y = FC), width = 0.2,         # Points
              color = "black", size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(unique(analysis_R_f$group_condition))
+
theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black")
      )+
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.1,  # space above highest point
           label = sig_label,
           size = 6))+ 
    scale_y_continuous(expand = c(0,0)) 

}
