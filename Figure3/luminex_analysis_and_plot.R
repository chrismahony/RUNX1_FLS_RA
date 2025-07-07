######### Read in data ############

library(readxl)
setwd("/rds/projects/c/croftap-runx1data01/luminex/repeat_jul7_2025")
# Read in and process data
Plate1_run1 <- read_excel("data_R.xlsx")

Plate1_run1 <- as.data.frame(Plate1_run1)

rownames(Plate1_run1) <- Plate1_run1$...1
Plate1_run1$...1 <- NULL

colnames(Plate1_run1) <- gsub("\\s+$", "", gsub("\\[.*?\\]|\\(.*?\\)", "", colnames(Plate1_run1)))

# Read in absolute standar conentration
standards <- read_excel("standards.xlsx")
colnames(standards) <- gsub("\\s+$", "", gsub("\\[.*?\\]|\\(.*?\\)", "", colnames(standards)))

standards <- as.data.frame(standards)
#rownames(standards) <- standards$...1
standards$...1 <- NULL



######### Process and intial plot ############


# Plot regression and calculate equation
standards_1_50 <- Plate1_run1 %>% filter(Condition == "Standard") %>% filter(Dilution == "1_50")

cols <- colnames(standards_1_50)[-c(1:3, 20:23)]

# Read in absolute standar conentration
standards <- read_excel("standards.xlsx")
colnames(standards) <- gsub("\\s+$", "", gsub("\\[.*?\\]|\\(.*?\\)", "", colnames(standards)))

standards <- as.data.frame(standards)
#rownames(standards) <- standards$...1
standards$...1 <- NULL



# Plot regression and calculate equation
standards_1_50 <- Plate1_run1 %>% filter(Condition == "Standard") %>% filter(Dilution == "1_50")

cols <- colnames(standards_1_50)[-c(1:3, 20:23)]

# Read in absolute standar conentration
standards <- read_excel("standards.xlsx")
colnames(standards) <- gsub("\\s+$", "", gsub("\\[.*?\\]|\\(.*?\\)", "", colnames(standards)))

standards <- as.data.frame(standards)
#rownames(standards) <- standards$...1
standards$...1 <- NULL



standards2_dil <- standards/50

models <- list()


expr <- Plate1_run1 %>% filter(Condition %in% c("Control", "RUNX1")) %>% filter(Dilution == "1_50")
library(drc)


for (i in cols) {
mediablank_mean <- Plate1_run1 %>% dplyr::select(i, Well, Condition, Dilution)
mediablank_mean <- mediablank_mean %>% filter(Condition == "media_blank" & Dilution == "1_50")
mediablank_mean <- mean(as.numeric(mediablank_mean[,1]))

standards_1_50$conc <- standards2_dil[[i]]

print(standards_1_50 %>% 
    ggplot(aes(y = standards_1_50[[i]], x = conc)) + 
    geom_point() + 
    geom_smooth(color = "blue", se=F) +  # Adding a best-fit line
    ggtitle(paste("Best-fit for", i)))


# Linear regression
#models[[i]] <- lm(standards_1_50[[i]] ~ conc, data = standards_1_50)
#intercept <- coef(models[[i]])[1]
#slope <- coef(models[[i]])[2]
#expr[[i]] <- as.numeric(expr[[i]])
#expr[[i]] <- expr[[i]] - mediablank_mean
#expr[[i]] <- (expr[[i]] - intercept)/slope


# 4/5PL curve
models[[i]] <- drm(standards_1_50[[i]] ~ conc, data = standards_1_50, fct = LL.5())
expr[[i]] <- as.numeric(expr[[i]])
expr[[i]] <- expr[[i]] - mediablank_mean

# Use the fitted 4PL or 5PL model to predict the concentration
expr[[i]] <- as.data.frame(expr[[i]])
expr[[i]] <- predict(models[[i]], newdata = expr[[i]])

plot(models[[i]], type = "all", main = "4PL Curve Fit", xlab = "Concentration", ylab = i, col = "blue")


print(ggplot(expr, aes(x=Condition, y=expr[[i]])) + 
  geom_boxplot(notch=F)+ggtitle(i))


}




df_avg <- expr %>%
  group_by(Condition, Dilution, Replicate) %>% 
  dplyr::summarise(across(everything(), mean, na.rm = TRUE))

df_avg <- df_avg %>%
  mutate(batch = case_when(
    Replicate %in% c(1, 2) ~ "batch1",
    Replicate %in% c(3, 4, 5) ~ "batch2",
    Replicate %in% c(6, 7) ~ "batch3",
    TRUE ~ "batch4"
  ))

plots_list <- list()

for (i in cols) {
  
  plots_list[[i]] <- ggplot(df_avg, aes(x = Condition, y = .data[[i]], fill = Condition)) +
    geom_boxplot(notch = FALSE) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_compare_means(method = "t.test", label = "p.signif", label.y = max(df_avg[[i]], na.rm = TRUE) * 1) +  # t-test with significance asterisks
    scale_fill_manual(values = c("Control" = "gray", "RUNX1" = "red")) +
    theme_ArchR() +
    ggtitle(i) +
    RotatedAxis() +
    NoLegend()

  print(plots_list[[i]])
}

######### Final plots ############

for (i in cols[-c(5,9)]) {
  
        df_avg_tmp <- df_avg[-c(6,19,20),]
  write.csv(df_avg_tmp, "/rds/projects/c/croftap-runx1data01/luminex/repeat_jul7_2025/luminex_source.csv")
        
    df_avg_tmp <- df_avg_tmp %>% dplyr::select(print(i))
    
    colnames(df_avg_tmp) <- c("Condition", "Dilution", "Values")
    
    summary_data <- df_avg_tmp %>%
  dplyr::group_by(Condition) %>%
  dplyr::summarise(
    Mean = mean(Values, na.rm = TRUE),
    SD = sd(Values, na.rm = TRUE)
  )
    
  ttest_res <- t.test(Values ~ Condition, data = df_avg_tmp)

  p_val <- ttest_res$p.value
sig_label <- case_when(
  p_val < 0.0001 ~ "****",
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max=max(df_avg_tmp$Values)

print(summary_data %>% 
ggplot( aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill=c("lightgrey", "red")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +        # Error bars
  geom_jitter(data = df_avg_tmp, aes(x = Condition, y = Values ), width = 0.2,         # Points
              color = "black", size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Bar Plot with Individual Points and SD",
       y = "Measured Value")+ggtitle(print(i))
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
           y = y_max + 0.01,  # space above highest point
           label = sig_label,
           size = 6))+ 
    scale_y_continuous(expand = c(0,0)     ) 



  }
