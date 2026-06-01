# ==============================================================================
# BACHELOR THESIS: Designing the Bridge: Transitional Government Composition
#                  and Post-Civil War Stability
# Author: Moritz-Tristan Mössner
# University of Zurich, Department of Political Science
# Supervisor: Dr. Belen Gonzalez
# ------------------------------------------------------------------------------
# SCRIPT 1 of 2: DATA MANAGEMENT
# Loads the raw datasets, filters and links them, codes the outcome, builds the
# inclusivity index and control variables, and saves the analysis-ready dataset.
# Run this script before 02_analysis.R.
# ==============================================================================

#set seed for replication
set.seed(125)

#prep libraries and datasets
library(tidyverse)
library(lubridate)
library(marginaleffects)
library(stargazer)
library(knitr)
library(kableExtra)

pax_path  <- "Datasets/pax_data_2257_agreements_v10.csv"
ucdp_path <- "Datasets/UcdpPrioConflict_v25_1.rds"

pax  <- read.csv(pax_path, stringsAsFactors = FALSE)
ucdp <- readRDS(ucdp_path)

dir.create("Exported_Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Exported_Tables", showWarnings = FALSE, recursive = TRUE)
#dataprep
max_ucdp_year <- max(ucdp$year)  # 2024

#filter for relevant agreements, link them, coding
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

#drop rightcensored cases
final_data <- results %>% filter(Outcome != "Censored")

#sanity check
print(table(final_data$Outcome))


#index building -> count subcomp, calculate score
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
                   ifelse(GRaSubs > 0, 1, 0)) / 3,
    
    Inclusivity_Index = (dim_military    * 1.5) +
      (dim_territorial * 1.5) +
      (dim_political   * 1.0) +
      (dim_economic    * 0.5) +
      (dim_civil       * 0.5),
    #restricted index
    Inclusivity_Restricted = ((dim_military * 1.5) + (dim_territorial * 1.5)) * (5/3)
  )

summary(final_data$Inclusivity_Index)
summary(final_data$Inclusivity_Restricted)

#prepare final dataset
final_data %>%
  select(dim_political, dim_military, dim_economic, dim_territorial, dim_civil) %>%
  summary() %>%
  print()


#centering index and equal weighted index
final_data <- final_data %>%
  mutate(
    Inclusivity_c            = Inclusivity_Index - mean(Inclusivity_Index, na.rm = TRUE),
    Inclusivity_Restricted_c = Inclusivity_Restricted - mean(Inclusivity_Restricted, na.rm = TRUE),
    Inclusivity_Equal        = dim_political + dim_military + dim_economic +
      dim_territorial + dim_civil,
    Inclusivity_Equal_c      = Inclusivity_Equal - mean(Inclusivity_Equal, na.rm = TRUE)
  )

cat("Cor(X, X^2) uncentered:",
    round(cor(final_data$Inclusivity_Index, final_data$Inclusivity_Index^2), 3), "\n")
cat("Cor(X, X^2) centered  :",
    round(cor(final_data$Inclusivity_c,    final_data$Inclusivity_c^2),    3), "\n")


#control variables
conflict_traits <- ucdp %>%
  arrange(conflict_id, year) %>%
  group_by(conflict_id) %>%
  summarise(
    start_date      = min(as.Date(start_date)),
    incompatibility = first(incompatibility),
    region          = first(region)
  )

#Conflict_Duration -> calendar years from onset to agreement
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
    Incompatibility_Label = factor(Incompatibility_Label,
                                   levels = c("Territorial",
                                              "Government",
                                              "Government & Territory"))
  )


#save analysis-ready datasets for 02_analysis.R
#final_data is the analysis sample; results is the pre-censoring object needed
#for the truncation check in the analysis script
saveRDS(final_data, "Datasets/final_data.rds")
saveRDS(results,    "Datasets/results.rds")
