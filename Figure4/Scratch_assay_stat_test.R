library(readxl)
measurments <- read_excel("/rds/projects/c/croftap-labdata2/Chris/data_plot_RUNX1_paper/figure4c.xlsx")


measurement_cols <- setdiff(
  names(measurments)[sapply(measurments, is.numeric)],
  c("Time", "Invasion")  # exclude grouping vars that are numeric
)

df_avg$Time <- paste("T", df_avg$Time, sep="")

df_avg <- measurments %>% 
  group_by(Time, Condition, Invasion) %>%
  dplyr::summarise(across(all_of(measurement_cols), mean, na.rm = TRUE), .groups = "drop")


res.aov2 <- aov(Invasion ~ Condition + Time, data = measurments)
summary(res.aov2)
plot(res.aov2, 1)
plot(res.aov2, 2)
aov_residuals <- residuals(object = res.aov2)
shapiro.test(x = aov_residuals )
TukeyHSD(res.aov2, which = "Condition")

library(purrr)
pairwise_results <- measurments %>%
  group_by(Time) %>%
  nest() %>%
  mutate(
    pw_tests = map(data, ~ {
      # Apply pairwise t-test
      result <- pairwise.t.test(.x$Invasion, .x$Condition, p.adjust.method = "none")
      
      # Tidy the result into a data frame
      broom::tidy(result)
    })
  ) %>%
  dplyr::select(Time, pw_tests) %>%
  unnest(pw_tests)
