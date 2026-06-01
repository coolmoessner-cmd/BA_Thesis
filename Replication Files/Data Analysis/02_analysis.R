# ==============================================================================
# BACHELOR THESIS: Designing the Bridge: Transitional Government Composition
#                  and Post-Civil War Stability
# Author: Moritz-Tristan MC6ssner
# University of Zurich, Department of Political Science
# Supervisor: Dr. Belen Gonzalez
# ------------------------------------------------------------------------------
# SCRIPT 2 of 2: ANALYSIS
# Loads the analysis-ready dataset produced by 01_data_management.R, runs the
# diagnostics, models and robustness checks, and exports all tables and figures.
# Run 01_data_management.R before this script.
# ==============================================================================

#set seed for replication
set.seed(125)

#prep libraries
library(tidyverse)
library(lubridate)
library(marginaleffects)
library(stargazer)
library(knitr)
library(kableExtra)

#ensure output folders exist
dir.create("Exported_Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Exported_Tables", showWarnings = FALSE, recursive = TRUE)

#load analysis-ready datasets from 01_data_management.R
final_data <- readRDS("Datasets/final_data.rds")
results    <- readRDS("Datasets/results.rds")


#separation diagnostic
cat("\n--- SEPARATION DIAGNOSTIC: Incompatibility x Failure ---\n")
print(table(
  Incompatibility = final_data$Incompatibility,
  Failed          = final_data$failed
))

#conflict duration by incompa
final_data %>%
  group_by(Incompatibility_Label) %>%
  summarise(
    n        = n(),
    mean_dur = round(mean(Conflict_Duration, na.rm = TRUE), 1),
    .groups  = "drop"
  ) %>%
  print()

#truncation check
failure_timing <- results %>%
  filter(Outcome == "Failure (Recurrence)") %>%
  rowwise() %>%
  mutate(
    first_failure_year = min(
      unlist(conflict_years)[
        unlist(conflict_years) >= window_start &
          unlist(conflict_years) <= window_end
      ],
      na.rm = TRUE
    ),
    years_to_failure = first_failure_year - Agre_Year
  ) %>%
  ungroup()

print(table(Years_to_Failure = failure_timing$years_to_failure))
cat("\nProportion failing in years 1-2:",
    round(mean(failure_timing$years_to_failure <= 2, na.rm = TRUE), 3), "\n")
cat("Proportion failing in years 4-5:",
    round(mean(failure_timing$years_to_failure >= 4, na.rm = TRUE), 3), "\n")


#lets do some modelling
model_logit <- glm(
  failed ~ Inclusivity_c + I(Inclusivity_c^2) +
    Conflict_Duration + Incompatibility,
  family = binomial(link = "logit"),
  data   = final_data
)

model_linear <- glm(
  failed ~ Inclusivity_c +
    Conflict_Duration + Incompatibility,
  family = binomial(link = "logit"),
  data   = final_data
)

lrt_result <- anova(model_linear, model_logit, test = "LRT")
print(lrt_result)

#territorial subsample
terr_data <- final_data %>% filter(incompatibility == 1)

model_terr <- glm(
  failed ~ Inclusivity_c + I(Inclusivity_c^2) + Conflict_Duration,
  family = binomial(link = "logit"),
  data   = terr_data
)

model_lpm <- lm(
  failed ~ Inclusivity_c + I(Inclusivity_c^2) +
    Conflict_Duration + Incompatibility,
  data = final_data
)

model_restricted <- glm(
  failed ~ Inclusivity_Restricted_c + I(Inclusivity_Restricted_c^2) +
    Conflict_Duration + Incompatibility,
  family = binomial(link = "logit"),
  data   = final_data
)

model_equal <- glm(
  failed ~ Inclusivity_Equal_c + I(Inclusivity_Equal_c^2) +
    Conflict_Duration + Incompatibility,
  family = binomial(link = "logit"),
  data   = final_data
)


#outputs of results
cat("\n--- PRIMARY MODEL ---\n")
summary(model_logit)

cat("\n--- LINEAR-ONLY MODEL (FOR LRT) ---\n")
summary(model_linear)

cat("\n--- ROBUSTNESS CHECK 1: TERRITORIAL SUBSAMPLE (n =", nrow(terr_data), ") ---\n")
summary(model_terr)

cat("\n--- ROBUSTNESS CHECK 2: LINEAR PROBABILITY MODEL ---\n")
summary(model_lpm)

cat("\n--- ROBUSTNESS CHECK 3: RESTRICTED INDEX ---\n")
summary(model_restricted)

cat("\n--- ROBUSTNESS CHECK 4: EQUAL-WEIGHTED INDEX ---\n")
summary(model_equal)

b1 <- coef(model_logit)["Inclusivity_c"]
b2 <- coef(model_logit)["I(Inclusivity_c^2)"]
mean_incl <- mean(final_data$Inclusivity_Index, na.rm = TRUE)

cat("\n--- FAILURE RATE BY INCLUSIVITY QUARTILE ---\n")
quartile_table <- final_data %>%
  mutate(Index_Quartile = ntile(Inclusivity_Index, 4)) %>%
  group_by(Index_Quartile) %>%
  summarise(
    n            = n(),
    mean_index   = round(mean(Inclusivity_Index), 2),
    failure_rate = round(mean(failed), 2),
    se           = round(sqrt(failure_rate * (1 - failure_rate) / n), 3)
  )
print(quartile_table)

cat("\n--- INDEX DISTRIBUTION BY BAND (Table 1 in thesis) ---\n")
band_table <- final_data %>%
  mutate(band = cut(Inclusivity_Index,
                    breaks = c(0, 1, 2, 3, 4, 5),
                    include.lowest = TRUE,
                    right = FALSE)) %>%
  group_by(band) %>%
  summarise(
    n            = n(),
    pct          = round(n() / nrow(final_data) * 100, 1),
    failure_rate = round(mean(failed) * 100, 1),
    min_index    = round(min(Inclusivity_Index), 2),
    max_index    = round(max(Inclusivity_Index), 2)
  )
print(band_table)


#Export Tables

#Table 2: main results
stargazer(model_logit,
          type  = "latex",
          title = "Main Results: Transitional Government Inclusivity and Conflict Recurrence (N = 44)",
          label = "tab:results_main",
          out   = "Exported_Tables/Table2_MainResults.tex",
          covariate.labels = c("Inclusivity Index (centred)",
                               "Inclusivity Index$^2$ (centred)",
                               "Conflict Duration",
                               "Incompatibility (Gov.)"),
          dep.var.labels = "Conflict Recurrence (binary)",
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "$^{*}$p $<$ .05; $^{**}$p $<$ .01; $^{***}$p $<$ .001. Two-tailed tests. Inclusivity terms mean-centred.",
          notes.append = FALSE)

#Table 3: robustness checks
stargazer(model_lpm, model_terr, model_restricted, model_equal,
          type  = "latex",
          title = "Robustness Checks: LPM, Territorial Subsample, Restricted Index, Equal Weights",
          label = "tab:results_robust",
          out   = "Exported_Tables/Table3_RobustnessChecks.tex",
          column.labels = c("LPM", "Terr. Sub.", "Restricted", "Equal Weights"),
          dep.var.labels = "Conflict Recurrence",
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "$^{*}$p $<$ .05; $^{**}$p $<$ .01; $^{***}$p $<$ .001. Two-tailed tests. All inclusivity terms mean-centred.",
          notes.append = FALSE)


#plots and exports

mean_incl_val <- mean(final_data$Inclusivity_Index)
sd_incl_val   <- sd(final_data$Inclusivity_Index)

#Figure 1: outcome bar chart
ggplot(final_data, aes(x = Outcome, fill = Outcome)) +
  geom_bar(color = "black", alpha = 0.8, width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Failure (Recurrence)" = "firebrick",
                               "Success (Peace)"      = "steelblue")) +
  labs(
    title    = "Empirical Record of Transitional Governments",
    subtitle = "Conflict Recurrence vs. Peace (5-Year Window)",
    x        = "Transition Outcome",
    y        = "Number of Interim Governments"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle   = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.text.x     = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 12)
  )
ggsave("Exported_Figures/Figure1_OutcomeBar.png",
       width = 8, height = 6, dpi = 300)


#Figure 4: regional distribution
ggplot(final_data, aes(x = Region, fill = Region)) +
  geom_bar(color = "black", alpha = 0.8, width = 0.6) +
  geom_text(stat = "count", aes(label = after_stat(count)),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_x_discrete(labels = c("1" = "Europe", "2" = "Middle East",
                              "3" = "Asia",   "4" = "Africa",
                              "5" = "Americas")) +
  labs(
    title    = "Regional Distribution of Transitional Governments",
    subtitle = "Intrastate Peace Agreements (1990-2019)",
    x        = "Geographic Region",
    y        = "Number of Agreements"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle   = element_text(hjust = 0.5, size = 12, color = "grey30"),
    axis.title      = element_text(size = 12)
  )
ggsave("Exported_Figures/Figure4_Regions.png",
       width = 8, height = 6, dpi = 300)


#Figure 5: year x inclusivity scatter
ggplot(final_data,
       aes(x = Agre_Year, y = Inclusivity_Index,
           color = Incompatibility_Label,
           size  = Conflict_Duration)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = mean_incl_val,
             linetype = "dashed", color = "grey30", linewidth = 0.6) +
  annotate("text", x = 2019, y = mean_incl_val + 0.12,
           label = "Mean", size = 3, color = "grey30", hjust = 1) +
  scale_color_manual(
    name   = "Incompatibility Type",
    values = c("Territorial"            = "steelblue",
               "Government"             = "firebrick",
               "Government & Territory" = "darkorange")
  ) +
  scale_size_continuous(
    name   = "Conflict Duration (years)",
    range  = c(2, 6),
    breaks = c(5, 10, 20, 30)
  ) +
  scale_x_continuous(breaks = seq(1990, 2019, 5)) +
  scale_y_continuous(limits = c(0, 5), breaks = 0:5) +
  labs(
    title    = "Transitional Governments by Year and Inclusivity Score",
    subtitle = paste0("Intrastate peace agreements 1990-2019 (N = ",
                      nrow(final_data), ")"),
    x        = "Year of Agreement",
    y        = "Weighted Inclusivity Index (0-5)"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title      = element_text(size = 11),
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9)
  )
ggsave("Exported_Figures/Figure5_Distribution.png",
       width = 9, height = 6, dpi = 300)


#Figure 6: predicted curve + rug
pred_data <- data.frame(
  Inclusivity_Index = seq(0, max(final_data$Inclusivity_Index),
                          length.out = 200),
  Conflict_Duration = mean(final_data$Conflict_Duration, na.rm = TRUE),
  Incompatibility   = factor(1, levels = levels(final_data$Incompatibility))
)
pred_data$Inclusivity_c <- pred_data$Inclusivity_Index -
  mean(final_data$Inclusivity_Index, na.rm = TRUE)
pred_data$predicted_prob <- predict(model_logit,
                                    newdata = pred_data,
                                    type    = "response")

ggplot() +
  annotate("rect",
           xmin = max(final_data$Inclusivity_Index) + 0.1, xmax = 5,
           ymin = 0, ymax = 1,
           alpha = 0.07, fill = "grey20") +
  annotate("text", x = 4.3, y = 0.5,
           label = "Empirically\nunobserved\nrange",
           size = 3, color = "grey40", fontface = "italic") +
  geom_line(data = pred_data,
            aes(x = Inclusivity_Index, y = predicted_prob),
            color = "darkred", linewidth = 1.2) +
  geom_rug(data = filter(final_data, failed == 1),
           aes(x = Inclusivity_Index),
           sides = "b", color = "firebrick",
           alpha = 0.7, length = unit(0.04, "npc")) +
  geom_rug(data = filter(final_data, failed == 0),
           aes(x = Inclusivity_Index),
           sides = "b", color = "steelblue",
           alpha = 0.7, length = unit(0.04, "npc")) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5)) +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Predicted Probability of Conflict Recurrence by Inclusivity",
    subtitle = paste0("Logistic regression, Duration and Incompatibility at mean/modal values (N = ",
                      nrow(final_data), ")"),
    x        = "Weighted Inclusivity Index (0-5)",
    y        = "Predicted Probability of Conflict Recurrence",
    caption  = "Red rug: failures. Blue rug: successes."
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title    = element_text(size = 11),
    plot.caption  = element_text(size = 9, color = "grey50")
  )
ggsave("Exported_Figures/Figure6_PredictedCurve.png",
       width = 9, height = 6, dpi = 300)


#Figure 7: substantive effect - incompatibility type
pred_incomp <- predictions(
  model_logit,
  newdata = datagrid(
    Incompatibility   = factor(c(1, 2),
                               levels = levels(final_data$Incompatibility)),
    Conflict_Duration = mean(final_data$Conflict_Duration, na.rm = TRUE),
    Inclusivity_c     = 0
  )
)

pred_incomp_df <- as.data.frame(pred_incomp) %>%
  mutate(
    Incompatibility_Label = ifelse(Incompatibility == 1,
                                   "Territorial", "Government")
  )

ggplot(pred_incomp_df,
       aes(x = Incompatibility_Label, y = estimate,
           ymin = conf.low, ymax = conf.high,
           color = Incompatibility_Label)) +
  geom_point(size = 4) +
  geom_errorbar(width = 0.12, linewidth = 0.9) +
  scale_color_manual(
    values = c("Territorial" = "steelblue", "Government" = "firebrick"),
    guide  = "none"
  ) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  labs(
    title    = "Predicted Recurrence Probability by Incompatibility Type",
    subtitle = "Inclusivity and Duration held at sample mean, with 95% confidence intervals",
    x        = "Incompatibility Type",
    y        = "Predicted Probability of Conflict Recurrence"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title    = element_text(size = 11)
  )
ggsave("Exported_Figures/Figure7_IncompEffect.png",
       width = 7, height = 6, dpi = 300)


#Figure 8: substantive effect - conflict duration
duration_range <- data.frame(
  Conflict_Duration = seq(min(final_data$Conflict_Duration, na.rm = TRUE),
                          max(final_data$Conflict_Duration, na.rm = TRUE),
                          length.out = 100),
  Incompatibility   = factor(1, levels = levels(final_data$Incompatibility)),
  Inclusivity_c     = 0
)

pred_duration    <- predictions(model_logit, newdata = duration_range)
pred_duration_df <- as.data.frame(pred_duration)

ggplot(pred_duration_df,
       aes(x = Conflict_Duration, y = estimate,
           ymin = conf.low, ymax = conf.high)) +
  geom_ribbon(alpha = 0.15, fill = "darkred") +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_rug(data = final_data,
           aes(x = Conflict_Duration),
           sides = "b", inherit.aes = FALSE,
           alpha = 0.4, color = "grey40") +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  labs(
    title    = "Substantive Effect of Conflict Duration on Recurrence Probability",
    subtitle = "Territorial incompatibility, inclusivity held at sample mean, with 95% CI",
    x        = "Pre-Agreement Conflict Duration (years)",
    y        = "Predicted Probability of Conflict Recurrence"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.title    = element_text(size = 11)
  )
ggsave("Exported_Figures/Figure8_DurationEffect.png",
       width = 9, height = 6, dpi = 300)


#case list export for appendix
case_list <- final_data %>%
  transmute(
    Country = Con,
    `Agreement Year` = Agre_Year,
    `UCDP Conflict ID` = UcdpID,
    `Incompatibility Type` = Incompatibility_Label,
    `Inclusivity Score` = round(Inclusivity_Index, 2),
    `Conflict Duration` = Conflict_Duration,
    `Recurrence within Five Years` = ifelse(failed == 1, "Yes", "No")
  ) %>%
  arrange(Country, `Agreement Year`)

case_list_latex <- case_list %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = "Case List of Transitional Governments Included in the Analysis",
    label = "tab:case_list"
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    font_size = 8
  )

save_kable(
  case_list_latex,
  file = "Exported_Tables/Appendix_Table_A1_CaseList.tex"
)

sessionInfo()
