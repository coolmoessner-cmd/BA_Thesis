---
title: "Replication: Designing the Bridge: Transitional Government Composition and Post-Civil War Stability"
subtitle: "Bachelor Thesis — University of Zurich, Department of Political Science"
author: "Moritz-Tristan Mössner"
date: "2026"
output:
  pdf_document:
    toc: true
    toc_depth: 3
    number_sections: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo    = TRUE,
  warning = FALSE,
  message = FALSE
)

# Set seed for replication
set.seed(125)
dir.create("Exported_Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Exported_Tables", showWarnings = FALSE, recursive = TRUE)
```

# Overview

This document replicates the empirical analysis reported in the BA thesis *Designing the Bridge: Transitional Government Composition and Post-Civil War Stability*. All code is self-contained. To reproduce the results, set the working directory to the `Replication Files` root folder and ensure both datasets are present in the `Datasets/` subfolder.

**Data sources:**

* PA-X Peace Agreement Database Version 10 (Bell & Badanjak, 2019)
* UCDP/PRIO Armed Conflict Dataset Version 25.1 (Davies et al., 2025; Gleditsch et al., 2002)

---

# Libraries and Data Loading

```{r libraries}
library(tidyverse)
library(lubridate)
library(marginaleffects)
library(stargazer)
library(knitr)
library(kableExtra)
```

The datasets are loaded using relative paths. Run this script from the `Replication Files` root folder.

```{r load-data}
pax_path <- "Datasets/pax_data_2257_agreements_v10.csv"
ucdp_path <- "Datasets/UcdpPrioConflict_v25_1.rds"

pax <- read.csv(pax_path, stringsAsFactors = FALSE)
ucdp <- readRDS(ucdp_path)

# Identify latest observable year for right-censoring cutoff
max_ucdp_year <- max(ucdp$year)
cat("Latest UCDP year:", max_ucdp_year, "\n")
```

---

# Data Preparation and Outcome Coding

PA-X agreements are filtered to intrastate transitional governments only (`Interim == "Yes"`, `Contp == "Government/territory"`). Each agreement is linked to its UCDP conflict ID. Conflict recurrence is then coded as a binary variable: 1 if armed conflict meeting the UCDP threshold recurs within a five-year observation window (`t+1` to `t+5`), 0 otherwise. Cases where the observation window extends beyond the UCDP data coverage are right-censored and dropped.

```{r outcome-coding}
results <- pax %>%
  filter(Interim == "Yes", Contp == "Government/territory") %>%
  mutate(
    Agre_Year = year(as.Date(Dat)),
    UcdpID    = as.numeric(str_extract(UcdpCon, "\\d+"))
  ) %>%
  filter(!is.na(UcdpID)) %>%
  rowwise() %>%
  mutate(
    conflict_years = list(ucdp$year[ucdp$conflict_id == UcdpID]),
    window_start   = Agre_Year + 1,
    window_end     = Agre_Year + 5,
    is_censored    = window_end > max_ucdp_year,
    recurrence     = any(unlist(conflict_years) >= window_start &
                           unlist(conflict_years) <= window_end),
    failed         = ifelse(recurrence, 1, 0),
    Outcome        = case_when(
      is_censored ~ "Censored",
      recurrence  ~ "Failure (Recurrence)",
      !recurrence ~ "Success (Peace)"
    )
  ) %>%
  ungroup()

# Drop right-censored cases
final_data <- results %>% filter(Outcome != "Censored")

cat("N after censoring filter:", nrow(final_data), "\n")
print(table(final_data$Outcome))
```

---

# Inclusivity Index Construction

## Dimension Scores

Each of the five theoretical dimensions is scored on a 0--1 scale by counting which PA-X subcomponents are present in the agreement and dividing by the total possible subcomponents in that category.

```{r index-dimensions}
final_data <- final_data %>%
  mutate(
    dim_political = (ifelse(PpsGe    > 0, 1, 0) +
                       ifelse(PpsEx    > 0, 1, 0) +
                       ifelse(PpsOro   > 0, 1, 0) +
                       ifelse(PpsOthPr > 0, 1, 0) +
                       ifelse(PpsVet   > 0, 1, 0) +
                       ifelse(PpsAut   > 0, 1, 0) +
                       ifelse(PpsInt   > 0, 1, 0) +
                       ifelse(PpsOth   > 0, 1, 0)) / 8,

    dim_military = (ifelse(MpsMe  > 0, 1, 0) +
                      ifelse(MpsJt  > 0, 1, 0) +
                      ifelse(MpsPro > 0, 1, 0) +
                      ifelse(MpsOth > 0, 1, 0)) / 4,

    dim_economic = (ifelse(EpsRes > 0, 1, 0) +
                      ifelse(EpsFis > 0, 1, 0) +
                      ifelse(EpsOth > 0, 1, 0)) / 3,

    dim_territorial = (ifelse(TpsSub > 0, 1, 0) +
                         ifelse(TpsLoc > 0, 1, 0) +
                         ifelse(TpsAut > 0, 1, 0) +
                         ifelse(TpsOth > 0, 1, 0)) / 4,

    dim_civil = (ifelse(Civso   > 0, 1, 0) +
                   ifelse(GeWom   > 0, 1, 0) +
                   ifelse(GRaSubs > 0, 1, 0)) / 3
  )

final_data %>%
  select(dim_political, dim_military, dim_economic, dim_territorial, dim_civil) %>%
  summary() %>%
  print()
```

## Weighted Composite Index

Dimensions are combined using theoretical weights. This produces a 0--5 scale. A restricted index using only the military and territorial dimensions is also constructed for robustness.

```{r index-composite}
final_data <- final_data %>%
  mutate(
    Inclusivity_Index = (dim_military    * 1.5) +
                        (dim_territorial * 1.5) +
                        (dim_political   * 1.0) +
                        (dim_economic    * 0.5) +
                        (dim_civil       * 0.5),

    Inclusivity_Restricted = ((dim_military * 1.5) +
                              (dim_territorial * 1.5)) * (5/3)
  )

cat("Full index distribution:\n")
summary(final_data$Inclusivity_Index)

cat("\nRestricted index distribution:\n")
summary(final_data$Inclusivity_Restricted)
```

## Mean-Centring and Equal-Weighted Index

To reduce collinearity between the linear and squared terms, the index is mean-centred before squaring (Aiken & West, 1991). An equal-weighted index, where all dimensions receive a weight of 1.0, is constructed for robustness checks.

```{r index-centering}
final_data <- final_data %>%
  mutate(
    Inclusivity_c            = Inclusivity_Index - mean(Inclusivity_Index, na.rm = TRUE),
    Inclusivity_Restricted_c = Inclusivity_Restricted - mean(Inclusivity_Restricted, na.rm = TRUE),
    Inclusivity_Equal        = dim_political + dim_military + dim_economic +
                               dim_territorial + dim_civil,
    Inclusivity_Equal_c      = Inclusivity_Equal - mean(Inclusivity_Equal, na.rm = TRUE)
  )

cat("Collinearity diagnostic:\n")
cat("Cor(X, X^2) uncentered:",
    round(cor(final_data$Inclusivity_Index, final_data$Inclusivity_Index^2), 3), "\n")
cat("Cor(X, X^2) centered  :",
    round(cor(final_data$Inclusivity_c, final_data$Inclusivity_c^2), 3), "\n")
```

---

# Control Variables

Conflict duration and incompatibility type are derived from the UCDP/PRIO Armed Conflict Dataset and put together into the analysis dataset.

```{r control-variables}
# Arrange before first() so the coding is deterministic per conflict_id
conflict_traits <- ucdp %>%
  arrange(conflict_id, year) %>%
  group_by(conflict_id) %>%
  summarise(
    start_date      = min(as.Date(start_date)),
    incompatibility = first(incompatibility),
    region          = first(region),
    .groups         = "drop"
  )

final_data <- final_data %>%
  left_join(conflict_traits, by = c("UcdpID" = "conflict_id")) %>%
  mutate(
    Conflict_Duration = Agre_Year - year(start_date),
    Incompatibility   = as.factor(incompatibility),
    Region            = as.factor(region),
    Incompatibility_Label = case_when(
      incompatibility == 1 ~ "Territorial",
      incompatibility == 2 ~ "Government",
      incompatibility == 3 ~ "Government & Territory"
    ),
    Incompatibility_Label = factor(
      Incompatibility_Label,
      levels = c("Territorial", "Government", "Government & Territory")
    )
  )
```

---

# Diagnostics

## Separation Diagnostic

The table below checks whether incompatibility type predicts failure perfectly in any category.

```{r separation}
cat("\n--- SEPARATION DIAGNOSTIC: Incompatibility x Failure ---\n")
print(table(
  Incompatibility = final_data$Incompatibility,
  Failed          = final_data$failed
))

cat("\nLabelled version:\n")
print(table(
  Incompatibility = final_data$Incompatibility_Label,
  Failed          = final_data$failed
))
```

```{r duration-by-incompatibility}
# Mean pre-agreement conflict duration by incompatibility type
# (supports the "government cases average ~34 years" statement in the thesis)
final_data %>%
  group_by(Incompatibility_Label) %>%
  summarise(
    n        = n(),
    mean_dur = round(mean(Conflict_Duration, na.rm = TRUE), 1),
    .groups  = "drop"
  ) %>%
  print()
```


## Truncation Check

This checks whether the five-year binary observation window introduces meaningful truncation bias by examining the timing of failures relative to the agreement date. 

```{r truncation}
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

cat("Years from agreement to first recurrence:\n")
print(table(Years_to_Failure = failure_timing$years_to_failure))
cat("\nProportion failing in years 1-2:",
    round(mean(failure_timing$years_to_failure <= 2, na.rm = TRUE), 3), "\n")
cat("Proportion failing in years 4-5:",
    round(mean(failure_timing$years_to_failure >= 4, na.rm = TRUE), 3), "\n")
```

---

# Regression Models

## Primary Model and LRT

The main model is a binary logistic regression with a quadratic polynomial term to test for the hypothesised U-shaped relationship.

```{r primary-model}
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

cat("Primary model results:\n")
summary(model_logit)

cat("\nLRT — quadratic vs. linear-only model:\n")
lrt_result <- anova(model_linear, model_logit, test = "LRT")
print(lrt_result)
```

## Robustness Checks

I estimate four robustness checks: (1) territorial subsample only, (2) linear probability model, (3) restricted military-territorial index, and (4) equal-weighted index.

```{r robustness}
# Territorial subsample removes the government incompatibility cases responsible for separation.
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

cat("\n--- ROBUSTNESS CHECK 1: TERRITORIAL SUBSAMPLE (n =", nrow(terr_data), ") ---\n")
summary(model_terr)

cat("\n--- ROBUSTNESS CHECK 2: LINEAR PROBABILITY MODEL ---\n")
summary(model_lpm)

cat("\n--- ROBUSTNESS CHECK 3: RESTRICTED INDEX ---\n")
summary(model_restricted)

cat("\n--- ROBUSTNESS CHECK 4: EQUAL-WEIGHTED INDEX ---\n")
summary(model_equal)
```

---

# Descriptive Tables

## Table 1 — Index Distribution by Band

```{r table1}
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
    max_index    = round(max(Inclusivity_Index), 2),
    .groups      = "drop"
  )

print(band_table)
```

## Failure Rate by Quartile

```{r quartile}
quartile_table <- final_data %>%
  mutate(Index_Quartile = ntile(Inclusivity_Index, 4)) %>%
  group_by(Index_Quartile) %>%
  summarise(
    n            = n(),
    mean_index   = round(mean(Inclusivity_Index), 2),
    failure_rate = round(mean(failed), 2),
    se           = round(sqrt(failure_rate * (1 - failure_rate) / n), 3),
    .groups      = "drop"
  )

print(quartile_table)
```

---

# Table Export

```{r export-tables}
dir.create("Exported_Tables", showWarnings = FALSE, recursive = TRUE)

# Table 2: main results
stargazer(
  model_logit,
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
  notes.append = FALSE
)

# Table 3: robustness checks
stargazer(
  model_lpm, model_terr, model_restricted, model_equal,
  type  = "latex",
  title = "Robustness Checks: LPM, Territorial Subsample, Restricted Index, Equal Weights",
  label = "tab:results_robust",
  out   = "Exported_Tables/Table3_RobustnessChecks.tex",
  column.labels = c("LPM", "Terr. Sub.", "Restricted", "Equal Weights"),
  dep.var.labels = "Conflict Recurrence",
  star.cutoffs = c(0.05, 0.01, 0.001),
  notes = "$^{*}$p $<$ .05; $^{**}$p $<$ .01; $^{***}$p $<$ .001. Two-tailed tests. All inclusivity terms mean-centred.",
  notes.append = FALSE
)
```

---

# Figures

```{r figure-setup}
dir.create("Exported_Figures", showWarnings = FALSE, recursive = TRUE)

mean_incl_val <- mean(final_data$Inclusivity_Index)
sd_incl_val   <- sd(final_data$Inclusivity_Index)
```

## Figure 1 — Empirical Record

```{r figure1, fig.width=8, fig.height=6}
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
```

## Figure 4 — Regional Distribution

```{r figure4-regions, fig.width=8, fig.height=6}
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
```

## Figure 5 — Distribution by Year and Inclusivity

```{r figure5-distribution, fig.width=9, fig.height=6}
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
```

## Figure 6 — Predicted Recurrence Curve

```{r figure6-predicted-curve, fig.width=9, fig.height=6}
pred_data <- data.frame(
  Inclusivity_Index = seq(0, max(final_data$Inclusivity_Index),
                          length.out = 200),
  Conflict_Duration = mean(final_data$Conflict_Duration, na.rm = TRUE),
  Incompatibility   = factor(1, levels = levels(final_data$Incompatibility))
)

pred_data$Inclusivity_c <- pred_data$Inclusivity_Index -
  mean(final_data$Inclusivity_Index, na.rm = TRUE)

pred_data$predicted_prob <- predict(
  model_logit,
  newdata = pred_data,
  type    = "response"
)

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
```

## Figure 7 — Incompatibility Type Effect

```{r figure7-incompatibility-effect, fig.width=7, fig.height=6}
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
```

## Figure 8 — Conflict Duration Effect

```{r figure8-duration-effect, fig.width=9, fig.height=6}
duration_range <- data.frame(
  Conflict_Duration = seq(min(final_data$Conflict_Duration, na.rm = TRUE),
                          max(final_data$Conflict_Duration, na.rm = TRUE),
                          length.out = 100),
  Incompatibility   = factor(1, levels = levels(final_data$Incompatibility)),
  Inclusivity_c     = 0
)

pred_duration <- predictions(model_logit, newdata = duration_range)
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
```

---

# Appendix Table A1 — Case List

The following code creates the case list used in the analysis and exports it as a LaTeX `longtable` for inclusion in the thesis appendix.

```{r appendix-case-list}
dir.create("Exported_Tables", showWarnings = FALSE, recursive = TRUE)

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

print(case_list)
```

---

# Session Information

```{r session-info}
sessionInfo()
```
