setwd("/rds/projects/c/croftap-stia-atac/CM_multiome/Functional_validation/mouse_tamoxifen_exp/exp2")
library(readxl)
scores <- read_excel("scores_avg.xlsx")

scores$D1 <- scores$D1 - scores$D0
scores$D2 <- scores$D2 - scores$D0
scores$D3 <- scores$D3 - scores$D0
scores$D4 <- scores$D4 - scores$D0
scores$D5 <- scores$D5 - scores$D0
scores$D6 <- scores$D6 - scores$D0
scores$D7 <- scores$D7 - scores$D0
scores$D8 <- scores$D8 - scores$D0

scores$D0 <- scores$D0 - scores$D0


#scores$D8 <- NULL

#scores <- scores %>% filter(Repeat == 2) 

# Include STIA scores from SAM (different Seurm batch)
#scores <- scores %>% filter(!mouse %in% c(2,6,7))

# Excluded STIA score from SAM
#scores <- scores %>% filter(mouse %in% c(1,3,4,5,8,9))


# Remove mouse where indicated was sick
scores <- scores %>% filter(!mouse %in% c(30))


measurement_cols <- grep("^D\\d+", names(scores), value = TRUE)

df_avg <- scores %>%filter(!treatment == "SAM") %>% 
  group_by(mouse, treatment, Joint) %>%
  dplyr::summarise(across(all_of(measurement_cols), mean, na.rm = TRUE), .groups = "drop")


#df_avg[df_avg < 0] <- 0


measurements2 <- df_avg  %>% pivot_longer(cols=colnames(df_avg)[-c(1:3)],
                    names_to='Day',
                    values_to='Score')


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(
      mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      n = sum(!is.na(x[[col]]))
    )
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Apply the updated function
df2 <- data_summary(measurements2, varname="Score", 
                    groupnames=c("treatment", "Day", "Joint"))

# Add SEM column
df2$sem <- df2$sd / sqrt(df2$n)



df2$Score_sd <- df2$Score + df2$sd


curve_plots <- list()
for (i in 1:length(unique(measurements2$Joint))) {
  curve_plots[[i]] <- df2 %>%
    filter(Joint == unique(measurements2$Joint)[[i]]) %>%
    ggplot(aes(x = Day, y = Score, group = treatment, color = treatment, shape = treatment)) +
    geom_line() +
    geom_errorbar(
      aes(ymax = Score + sem, ymin = Score),
      width = 0.3,
      position = position_dodge(0.01),
      color = "#4D4D4D"
    ) +
    geom_point(size = 4.5) +  # shape is inherited from aes()
    labs(title = "Runx1 mouse exp1", x = "STIA Day", y = "Measurement") +
    theme_classic() +
    scale_color_manual(values = c("CTRL" = '#4D4D4D', "TAM" = 'red', "SAM" = 'blue')) +
    scale_shape_manual(values = c("CTRL" = 16, "TAM" = 15, "SAM" = 15)) +  # circle for control, square for others
    ggtitle(unique(measurements2$Joint)[[i]])

  print(curve_plots[[i]])
}
measurements2$treatmant_day <- paste(measurements2$treatment, measurements2$Day, sep="_")



#curve_plots <- list()
#for (i in 1:length(unique(scores$Joint))){
#curve_plots[[i]] <- scores %>% filter(Joint == unique(scores$Joint)[[i]]) %>% 
#ggplot( aes(x=Day, y=Score, group=treatment, color=treatment)) + 
#  geom_line() +
#  geom_point()+
#  geom_errorbar(aes(ymax=Score_sd, ymin=Score), width=.2,
#                 position=position_dodge(0.05))+labs(title="Runx1 mouse exp1", x="STIA Day", y = "Measurement")+
#   theme_classic() +
#   scale_color_manual(values=c('#4D4D4D','red', 'blue'))+ggtitle(unique(scores$Joint)[[i]])#
#
#print(curve_plots[[i]])
#}



R_paw_AUC <- list()
for ( i in 1:length(unique(measurements2$mouse))){

measurements3 <- measurements2 %>% filter(Joint == "R_paw")

measurements4 <- measurements3 %>% filter(mouse == unique(measurements3$mouse)[[i]])

library(pracma)
x_numeric <- as.numeric(gsub("D", "", measurements4$Day))
R_paw_AUC[[i]] <- trapz(x_numeric, measurements4$Score)

}

R_paw_AUC <- R_paw_AUC %>% as.data.frame() %>% t %>% as.data.frame()
rownames(R_paw_AUC) <- unique(measurements2$mouse)

R_paw_AUC$mouse <- rownames(R_paw_AUC)

index <- match(R_paw_AUC$mouse, measurements3$mouse)
R_paw_AUC$treatment <- measurements3$treatment[index]


# Stats
measurements3 <- measurements2 %>% filter(Joint == "ankle")
res.aov2 <- aov(Score ~ treatment + Day, data = measurements3)
summary(res.aov2)
plot(res.aov2, 1)
plot(res.aov2, 2)
aov_residuals <- residuals(object = res.aov2)
shapiro.test(x = aov_residuals )
TukeyHSD(res.aov2, which = "treatment")

library(purrr)
pairwise_results <- measurements3 %>%
  group_by(Day) %>%
  nest() %>%
  mutate(
    pw_tests = map(data, ~ {
      # Apply pairwise t-test
      result <- pairwise.t.test(.x$Score, .x$treatment, p.adjust.method = "none")
      
      # Tidy the result into a data frame
      broom::tidy(result)
    })
  ) %>%
  dplyr::select(Day, pw_tests) %>%
  unnest(pw_tests)




measurements3 <- measurements2 %>% filter(Joint == "R_paw")
res.aov2 <- aov(Score ~ treatment + Day, data = measurements3)
summary(res.aov2)
plot(res.aov2, 1)
plot(res.aov2, 2)
aov_residuals <- residuals(object = res.aov2)
shapiro.test(x = aov_residuals )
TukeyHSD(res.aov2, which = "treatment")

library(purrr)
pairwise_results <- measurements3 %>%
  group_by(Day) %>%
  nest() %>%
  mutate(
    pw_tests = map(data, ~ {
      # Apply pairwise t-test
      result <- pairwise.t.test(.x$Score, .x$treatment, p.adjust.method = "none")
      
      # Tidy the result into a data frame
      broom::tidy(result)
    })
  ) %>%
  dplyr::select(Day, pw_tests) %>%
  unnest(pw_tests)


library(ggpubr)
library(rstatix)

stat.test <- R_paw_AUC  %>%
  t_test(V1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "treatment")

p1 <- ggplot(R_paw_AUC, aes(x=treatment, y=V1)) +
  geom_boxplot(fill=c('darkgrey',"darkred"))+geom_point()+
  theme_classic()+ stat_pvalue_manual(stat.test)+labs(title="R paw", x="", y = "AU")

p1





Ankle_AUC <- list()
for ( i in 1:length(unique(measurements2$mouse))){

measurements3 <- measurements2 %>% filter(Joint == "ankle")

measurements4 <- measurements3 %>% filter(mouse == unique(measurements3$mouse)[[i]])

library(pracma)
x_numeric <- as.numeric(gsub("D", "", measurements4$Day))
Ankle_AUC[[i]] <- trapz(x_numeric, measurements4$Score)

}

Ankle_AUC <- Ankle_AUC %>% as.data.frame() %>% t %>% as.data.frame()
rownames(Ankle_AUC) <- unique(measurements2$mouse)

Ankle_AUC$mouse <- rownames(Ankle_AUC)

index <- match(Ankle_AUC$mouse, measurements3$mouse)
Ankle_AUC$treatment <- measurements3$treatment[index]


library(ggpubr)
library(rstatix)

#Ankle_AUC <- Ankle_AUC %>% filter(mouse != c(11))

stat.test <- Ankle_AUC  %>%
  t_test(V1 ~ treatment) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "treatment")

p2 <- ggplot(Ankle_AUC, aes(x=treatment, y=V1)) +
  geom_boxplot(fill=c('darkgrey',"darkred"))+geom_point()+
  theme_classic()+ stat_pvalue_manual(stat.test)+labs(title="Ankle", x="", y = "AU")

p2
