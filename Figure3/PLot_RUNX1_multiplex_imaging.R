RUNX1_intenst$RUNX1 <- as.double(RUNX1_intenst$RUNX1)
 
RUNX1_intenst %>% filter(cat %in% c("low", "high")) %>% 
ggplot(aes(x=cat, y=RUNX1)) +stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
        geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")+geom_point(size=2.5)+theme_ArchR()



df <- RUNX1_intenst %>% filter(cat %in% c("low", "high"))
t.test(RUNX1 ~ cat, data=df)


summary_data <- RUNX1_intenst %>%
  dplyr::group_by(cat) %>%
  dplyr::summarise(
    Mean = mean(RUNX1),
    SD = sd(RUNX1)
  )

RUNX1_intenst <- RUNX1_intenst %>% filter(cat!="LL")
RUNX1_intenst <- RUNX1_intenst %>% filter(cat!="Aggr control")


t_res <- t.test(RUNX1 ~ cat, data = RUNX1_intenst)

p_val <- t_res$p.value
sig_label <- case_when(
  p_val < 0.001 ~ "***",
  p_val < 0.01  ~ "**",
  p_val < 0.05  ~ "*",
  TRUE          ~ "ns"
)

y_max <- max(RUNX1_intenst$RUNX1, na.rm = TRUE)


summary_data$cat <- factor(summary_data$cat, levels = c("low", "high"))
RUNX1_intenst$Condition <- factor(RUNX1_intenst$cat, levels = c("Unstimulated", "DLL4 treated"))

summary_data %>%  
ggplot( aes(x = cat, y = Mean)) +
  geom_bar(stat = "identity", color = "black", width = 0.8, fill=c("red", "lightgrey")) +               # Bars
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, size = 0.8) +        # Error bars
  geom_jitter(data = RUNX1_intenst, aes(x = cat, y = RUNX1), width = 0.15,         # Points
              color = "black", size = 3.5, alpha = 0.8) +
   theme(
      panel.background = element_blank(),       # remove inner panel background
      plot.background = element_blank(),        # remove outer background
      panel.grid.major = element_blank(),       # remove major grid lines
      panel.grid.minor = element_blank(),       # remove minor grid lines
      axis.line = element_line(color = "black"), # keep axis lines
    plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
  annotate("text", 
           x = 1.5,  # midpoint between bar 1 and 2
           y = y_max + 0.5,  # space above highest point
           label = sig_label,
           size = 6)
