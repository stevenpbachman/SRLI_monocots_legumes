
# Threat proportion by group ----

# first table is collated categories per group
# SRLI original (2010)
# SRLI current (2000-2010-2020)
# Red List (2010-2020) 2010 has 12K
# Predictions (2024)

# second table is collated % threatened per group

# save tables as csv and combine as excel sup mat later

#TODO add n to chart, but tricky with labels?
#TODO add RL 2024 to align with predictions?

library(rlang)
library(glue)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(writexl)
library(readxl)
library(scales)
library(patchwork)

work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir,"00_raw_data/SRLI_2026.csv"))
preds <- read.csv(paste0(work_dir,"00_raw_data/Table_S7_final_predictions.csv"))
RLI_stat <- SRLI_2026 %>% filter(SRLI_2024_analysis == "include") # only the included species
RLI_stat_thin <- RLI_stat %>% select(SRLI_2000_final_cat,
                                     SRLI_2010_final_cat,
                                     SRLI_2020_final_cat,
                                     WCVP_current_id,
                                     group,
                                     scientificName)

# Functions ----
# Red List category mapping
iucn_map <- c(
  "Critically Endangered" = "CR",
  "Data Deficient" = "DD",
  "Endangered" = "EN",
  "Extinct" = "EX",
  "Extinct in the Wild" = "EW",
  "Least Concern" = "LC",
  "Lower Risk/conservation dependent" = "LR/cd",
  "Lower Risk/least concern" = "LR/lc",
  "Lower Risk/near threatened" = "LR/nt",
  "Near Threatened" = "NT",
  "Vulnerable" = "VU"
)

#' Generate overall threat level using DD uncertainty
#'
#' @param categories dataframe with long format categories e.g. one row per species
#' @param count_field if you already have summarised data, define the field with the total per category
#'
#' @return data frame with lower, mid-point "best" and upper estimates
#' @export
#'
#' @examples
RL_threat_est <- function(categories, cat_field, count_field = NULL) {

  # Convert cat_field to symbol
  cat_sym <- ensym(cat_field)

  # standardise legacy IUCN categories
  categories <- categories %>%
    mutate(
      !!cat_sym := toupper(!!cat_sym),  # optional but safer
      !!cat_sym := recode(
        !!cat_sym,
        "LR/CD" = "NT",
        "LR/NT" = "NT",
        "LR/LC" = "LC",
        .default = !!cat_sym
      )
    )

  # Define expected categories as a dataframe
  all_categories <- data.frame(category = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD", "NE"))

  # Check for unexpected categories
  input_categories <- unique(categories[[as_string(cat_sym)]])
  unexpected <- setdiff(input_categories, all_categories$category)
  if (length(unexpected) > 0) {
    warning("Unexpected categories detected: ", paste(unexpected, collapse = ", "))
  }

  # Handle different input formats
  if (is.null(count_field)) {
    # Long form data - count occurrences of each category
    category_counts <- categories %>%
      filter(!is.na(!!cat_sym)) %>% #& !!cat_sym != "NE") %>%
      count(!!cat_sym, name = "n") %>%
      rename(category = !!cat_sym)
  } else {
    # Summary format - use provided count field
    count_sym <- ensym(count_field)
    category_counts <- categories %>%
      filter(!is.na(!!cat_sym)) %>% #& !!cat_sym != "NE") %>%
      filter(!is.na(!!count_sym)) %>%
      select(category = !!cat_sym, n = !!count_sym)
  }

  # Calculate total assessed (excluding EX)
  assessed <- category_counts %>%
    filter(!!cat_sym != "NE") %>%
    summarise(total = sum(n, na.rm = TRUE)) %>%
    pull(total)

  # Create complete category summary with zeros for missing categories
  categories_summary <- all_categories %>%
    left_join(category_counts, by = "category") %>%
    mutate(value = ifelse(is.na(n), 0, n)) %>%
    select(category, value)

  # Extract individual category counts
  EX <- categories_summary %>% filter(category == "EX") %>% pull(value)
  EW <- categories_summary %>% filter(category == "EW") %>% pull(value)
  CR <- categories_summary %>% filter(category == "CR") %>% pull(value)
  EN <- categories_summary %>% filter(category == "EN") %>% pull(value)
  VU <- categories_summary %>% filter(category == "VU") %>% pull(value)
  NT <- categories_summary %>% filter(category == "NT") %>% pull(value)
  LC <- categories_summary %>% filter(category == "LC") %>% pull(value)
  DD <- categories_summary %>% filter(category == "DD") %>% pull(value)
  NE <- categories_summary %>% filter(category == "NE") %>% pull(value)

  # Calculate threat estimates
  assessed_eval_only <- assessed        # standard IUCN
  assessed_all_known <- assessed + NE   # NE treated as part of total pool

  threatened <- EW + CR + EN + VU

  unknown_DD <- DD
  unknown_NE <- DD + NE

  # STANDARD IUCN (evaluated only)
  lower_bound <- threatened / (assessed_eval_only - EX)
  mid_point <- threatened /
    ((assessed_eval_only - EX) - unknown_DD)
  upper_bound <- (threatened + unknown_DD) /
    (assessed_eval_only - EX)

  # NE INCLUDED (full species pool)
  NE_lower <- threatened / (assessed_all_known - EX)
  NE_mid <- threatened /
    ((assessed_all_known - EX) - unknown_NE)
  NE_upper <- (threatened + unknown_NE) /
    (assessed_all_known - EX)


  # lower_bound <- (EW + CR + EN + VU) / (assessed - EX)
  # mid_point <- (EW + CR + EN + VU) / ((assessed - EX) - DD)
  # upper_bound <- (EW + CR + EN + VU + DD) / (assessed - EX)
  #
  # # NE INCLUDED (DD + NE unknown)
  # unknown_NE = DD + NE
  # NE_lower = (EW + CR + EN + VU) / (assessed - EX)
  # NE_mid = (EW + CR + EN + VU) / ((assessed - EX) - unknown_NE)
  # NE_upper = (EW + CR + EN + VU + unknown_NE) / (assessed - EX)

  percent_threatened <- data.frame(
    lower_bound = lower_bound,
    mid_point = mid_point,
    upper_bound = upper_bound
    #NE_lower = NE_lower,
    #NE_mid_point = NE_mid,
    #NE_upper

  )

  return(percent_threatened)
}

#' Calculate Red List Index value
#'
#' @param data_sufficient_counts
#' @param boot
#' @param year
#' @param group
#' @param dd_species
#'
#' @return df with RLI value and year, group
#' @export
#'
#' @examples
# Function to calculate RLI
rli_boot <- function(data_sufficient_counts, boot = FALSE, year, group, dd_species) {
  weights <- c(LC = 0, NT = 1, VU = 2, EN = 3, CR = 4, EW = 4, EX = 5)
  filtered_counts <- data_sufficient_counts[names(data_sufficient_counts) %in% names(weights)]
  #aligned_weights <- weights[names(aligned_counts)]
  #total_weight <- sum(aligned_counts * aligned_weights)
  #filtered_counts <- data_sufficient_counts[setdiff(names(data_sufficient_counts), c("DD", "NE"))]
  aligned_weights <- weights[names(filtered_counts)]
  total_weight <- sum(filtered_counts * aligned_weights)
  #total_weight <- sum(data_sufficient_counts * weights[names(data_sufficient_counts)])
  max_weight <- sum(filtered_counts) * 5
  result <- 1 - (total_weight / max_weight)

  if (boot == TRUE) {
    # run the bootstrapping
    # Bootstrapping parameters
    n_bootstrap <- 1000
    bootstrap_rli <- numeric(n_bootstrap)
    set.seed(123)  # For reproducibility
    category_proportions <- proportions(filtered_counts)

    # Bootstrapping loop
    for (i in 1:n_bootstrap) {
      # Randomly assign DD species to categories based on proportions
      dd_counts <- rmultinom(1, size = dd_species, prob = category_proportions)
      #names(dd_counts) <- names(weights)  # Assign names to match category_weights

      # Combine DD counts with data-sufficient counts
      total_counts <- filtered_counts + as.numeric(dd_counts)

      # Calculate RLI for this bootstrap iteration
      bootstrap_rli[i] <- rli_boot(total_counts, boot = FALSE, year, group, dd_species)
      # bs_total_weight <- sum(total_counts * weights)
      # bs_max_weight <- sum(data_sufficient_counts) * 5
      # bootstrap_rli[i] <- 1 - (bs_total_weight / bs_max_weight)
    }

    # Calculate confidence intervals
    ci_lower <- quantile(bootstrap_rli, probs = 0.025)
    ci_upper <- quantile(bootstrap_rli, probs = 0.975)

    result <- data.frame(
      Year = year,
      Group = group,
      RLI = result,
      ci_lower = ci_lower,
      ci_upper = ci_upper
    )
  }

  return(result)
}


# START HERE TO MAKE FIGURES ----
# Table S1. categories per group ####
# 1. SRLI current version ----
# gives the counts per category and by year
# make data long first
RLI_year_cats <- RLI_stat_thin %>%
  pivot_longer(
    cols = c(SRLI_2000_final_cat, SRLI_2010_final_cat, SRLI_2020_final_cat),
    names_to = "year",
    values_to = "category"
  ) %>%
  mutate(year = str_extract(year, "\\d{4}")) # Extract just the year

# then summarise by group and year
SRLI_2020_categories <- RLI_year_cats %>%
  group_by(year, group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "SRLI 2020")

# summarised as percentage threatened using the IUCN approach
SRLI_2020_perc_threat <- SRLI_2020_categories %>%
  group_by(year, group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "SRLI 2020")

# 1.b SRLI 2020 with preds ----
# counts are summarised from the SRLI 2026 excel
RLI_year_cats <- RLI_stat %>%
  pivot_longer(
    cols = c(SRLI_2000_final_cat, SRLI_2010_final_cat, SRLI_2020_final_cat),
    names_to = "year",
    values_to = "category"
  ) %>%
  mutate(year = str_extract(year, "\\d{4}")) # Extract just the year

RLI_year_cats <- RLI_year_cats %>%
  mutate(
    category = case_when(
      category == "DD" & predicted_threat == "threatened"     ~ "CR",
      category == "DD" & predicted_threat == "not threatened" ~ "LC",
      category == "NE" & predicted_threat == "threatened"     ~ "CR",
      category == "NE" & predicted_threat == "not threatened" ~ "LC",
      category == "NE" & is.na(predicted_threat) ~ "DD",
      is.na(category) & is.na(predicted_threat)              ~ "DD",
      TRUE ~ category
    )
  )

# then summarise by group and year
SRLI_2020_categories_with_preds <- RLI_year_cats %>%
  group_by(year, group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "SRLI 2020")

SRLI_2020_perc_threat_augmented <- SRLI_2020_categories_with_preds %>%
  group_by(year, group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "SRLI 2020 +preds")

# 2. SRLI original version ----
srli_orig_mon <- read_excel(paste0(work_dir,"00_raw_data/CombinedSpeciesList_SRLI_GreenPlantsRed_Paper.xlsx"),
                        sheet = "Monocots") %>%
  dplyr::mutate(group = "Monocots")

srli_orig_leg <- read_excel(paste0(work_dir,"00_raw_data/CombinedSpeciesList_SRLI_GreenPlantsRed_Paper.xlsx"),
                            sheet = "Legumes") %>%
  dplyr::mutate(group = "Legumes")

# combine
srli_orig <- bind_rows(srli_orig_mon, srli_orig_leg)

# now get the categories
SRLI_2010_categories <- srli_orig %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "SRLI 2010", year = "2010")

# get the %
SRLI_2010_perc_threat <- SRLI_2010_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "SRLI 2010", year = "2010")

# 3. Red List ----
# include 2010 and 2020
# 2010 ----
RL2010_3 <- read_excel(paste0(work_dir,"00_raw_data/plants on red list 2010_3.xls"))
RL2010_3_legs <- RL2010_3 %>% dplyr::filter(Family == "LEGUMINOSAE") %>%
  mutate(group = "Legumes")
RL2010_3_mons <- RL2010_3 %>% dplyr::filter(Phylum == "TRACHEOPHYTA",
                                            Class == "LILIOPSIDA") %>%
  mutate(group = "Monocots")

# combine
RL2010_3_comb <- bind_rows(RL2010_3_legs, RL2010_3_mons)
RL2010_3_comb$category <- RL2010_3_comb$`Red List status`

# now get the categories
RL_2010_categories <- RL2010_3_comb %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "Red List", year = "2010")

#RL_2010_categories <- bind_rows(RL_2010_categories, mons, legs)

# get the %
RL_2010_perc_threat <- RL_2010_categories %>%
  group_by(group) %>%
  group_map(~ {
    #assessed <- dplyr::filter(.x, !category %in% "NE")
    assessed <- dplyr::filter(.x)
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List", year = "2010")

#2015 ----
RL_2015 <- read.csv(paste0(work_dir,"00_raw_data/rl_assessment_details_20201215.csv"))
RL_2015 <- RL_2015 %>% filter(published_year <= 2015)

RL2015_legs <- RL_2015 %>% dplyr::filter(family == "FABACEAE") %>%
  mutate(group = "Legumes")
RL2015_mons <- RL_2015 %>% dplyr::filter(phylum == "TRACHEOPHYTA",
                                         class == "LILIOPSIDA") %>%
  mutate(group = "Monocots")

# combine
RL2015_comb <- bind_rows(RL2015_legs, RL2015_mons)

# now get the catgories
RL_2015_categories <- RL2015_comb %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "Red List", year = "2015")

# get the %
RL_2015_perc_threat <- RL_2015_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List", year = "2015")

# 2020 ----
RL_2020 <- read_excel(paste0(work_dir,"00_raw_data/RL_2020.xlsx"))

RL2020_legs <- RL_2020 %>% dplyr::filter(family == "FABACEAE") %>%
  mutate(group = "Legumes")
RL2020_mons <- RL_2020 %>% dplyr::filter(phylum == "TRACHEOPHYTA",
                                            class == "LILIOPSIDA") %>%
  mutate(group = "Monocots")

# combine
RL2020_comb <- bind_rows(RL2020_legs, RL2020_mons)

# now get the threat cats
RL_2020_categories <- RL2020_comb %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "Red List", year = "2020")

# get the %
RL_2020_perc_threat <- RL_2020_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List", year = "2020")


# 2024 ----
RL_2024 <- read.csv(paste0(work_dir,"00_raw_data/WCVPv13-RL2024-2.csv"))

RL2024_legs <- RL_2024 %>% dplyr::filter(familyName == "FABACEAE") %>%
  mutate(group = "Legumes")
RL2024_mons <- RL_2024 %>% dplyr::filter(phylumName == "TRACHEOPHYTA",
                                         className == "LILIOPSIDA") %>%
  mutate(group = "Monocots")

# combine
RL2024_comb <- bind_rows(RL2024_legs, RL2024_mons)

RL2024_comb <- RL2024_comb %>%
  mutate(category = recode(redlistCategory, !!!iucn_map))

# now get the threat cats
RL_2024_categories <- RL2024_comb %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "Red List", year = "2024")

# get the %
RL_2024_perc_threat <- RL_2024_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List", year = "2024")

# 2025 ----
RL_2025 <- read.csv(paste0(work_dir,"00_raw_data/WCVPv14-RL2025-2.csv"))

RL2025_legs <- RL_2025 %>% dplyr::filter(familyName == "FABACEAE") %>%
  mutate(group = "Legumes")
RL2025_mons <- RL_2025 %>% dplyr::filter(phylumName == "TRACHEOPHYTA",
                                         className == "LILIOPSIDA") %>%
  mutate(group = "Monocots")

# combine
RL2025_comb <- bind_rows(RL2025_legs, RL2025_mons)

RL2025_comb <- RL2025_comb %>%
  mutate(category = recode(redlistCategory, !!!iucn_map))

# now get the threat cats
RL_2025_categories <- RL2025_comb %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "Red List", year = "2025")

# get the %
RL_2025_perc_threat <- RL_2025_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List", year = "2025")

# 3b. Red List + preds ----
# merge the 2024 published red list with predictions
RL2024_comb <- RL2024_comb %>%
  dplyr::select(group, category, ipni_id) # red list mons + legs with correct categories

# add the groups
preds <- preds %>%
  mutate(group = case_when(
    family == "Fabaceae" ~ "Legumes",
  ))

leg_count <- preds %>%
  filter(!is.na(group)) %>%
  count()

# mons
APGIV_families <- read.csv(paste0(work_dir,"00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))

preds <- preds %>%
  left_join(APGIV_families, by = c(family = 'APG.IV.Family'))

preds <- preds %>%
  mutate(group = if_else(HigherGroups == "Monocots", "Monocots", group))

preds %>%
  count(group)

preds_mon_legs <- preds %>% dplyr::filter(group %in% c("Monocots", "Legumes"))

preds_mon_legs <- preds_mon_legs %>%
  mutate(
    category = case_when(
      is.na(category) & predicted_threat == "threatened"     ~ "CR",
      is.na(category) & predicted_threat == "not threatened" ~ "LC",
      is.na(category) & is.na(predicted_threat)              ~ "DD",
      TRUE ~ category
    )
  )

# now get the threat cats
RL_2024_preds_categories <- preds_mon_legs %>%
  group_by(group, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(complete.cases(category )) %>%
  mutate(dataset = "RL_2024+preds", year = "2024")

# get the %
RL_2024_preds_perc_threat <- RL_2024_preds_categories %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List 2024 +preds", year = "2024")

# 3.c Red List with DD ----
RL_2024_categories_with_DD <- RL_2024_categories %>%
  mutate(count = case_when(
    group == "Monocots" & category == "DD" ~ count + 68204,
    group == "Legumes" & category == "DD" ~ count + 15899,
    TRUE ~ count
  )) %>%
  mutate(dataset = "Red List 2024 +DD")

# get the %
RL_2024_DD_perc_threat <- RL_2024_categories_with_DD %>%
  group_by(group) %>%
  group_map(~ {
    assessed <- dplyr::filter(.x, !category %in% "NE")
    bind_cols(
      .y,
      tibble(n = sum(assessed$count)),
      RL_threat_est(.x, "category", "count")
    )
  }) %>%
  bind_rows() %>%
  mutate(dataset = "Red List 2024 +NE", year = "2024")


###### update this bit with the same red list DD method, not the preds paper
# 4. predictions ----
bart_summary <- read_excel(paste0(work_dir,"00_raw_data/bart-summary.xlsx"),
                           sheet = "pred-family")

pred_legs <- bart_summary %>%
  filter(family == "Fabaceae")

pred_legs <- pred_legs %>%
  rename(group = family,
         mid_point = threatened,
         lower_bound = .lower,
         upper_bound = .upper
           ) %>%
  select(-.width, -.interval, -.point)
pred_legs$year <- "2024"
pred_legs$dataset <- "prediction"
pred_legs$group <- "Legumes"
#pred_legs$n <- 22500

##predictors_angiosperm <- read.csv("C:/Users/sb42kg/Downloads/cluster_results/predictors-angiosperm-20230604-105231.csv")

# map monocot families to preds data
APGIV_families <- read.csv(paste0(work_dir,"00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))

pred_mono <- bart_summary %>%
  left_join(APGIV_families, by = c(family = 'APG.IV.Family')) %>%
  filter(HigherGroups == 'Monocots') %>%
  select(family, threatened, .lower, .upper)

pred_mons <- data.frame(year = "2024",
                       group = "Monocots",
                       lower_bound = mean(pred_mono$.lower),
                       mid_point = mean(pred_mono$threatened),
                      upper_bound = mean(pred_mono$.upper),
                      dataset = "prediction")

Preds_2024_perc_threat <- bind_rows(pred_legs, pred_mons)

# preds new method ----
# legumes
preds_raw <- read.csv(paste0(work_dir,"00_raw_data/Table_S7_final_predictions.csv"))

legs_preds <- preds_raw %>% filter(family == "Fabaceae", complete.cases(combined))
legs_summary <- legs_preds %>% group_by(combined, confidence) %>% count()

N <- sum(legs_summary$n)

# Extract counts
threat_conf  <- legs_summary$n[legs_summary$combined == "threatened" &
                                 legs_summary$confidence == "confident"]

threat_low<- legs_summary$n[legs_summary$combined == "threatened" &
                             legs_summary$confidence == "low_confidence"]

not_conf     <- legs_summary$n[legs_summary$combined == "not threatened" &
                             legs_summary$confidence == "confident"]

# Calculations
lower_bound  <- threat_conf / N
upper_bound  <- (threat_conf + threat_low) / N
best_est     <- threat_conf / (threat_conf + not_conf)

# Combine results
RL_2024_preds_perc_threat_legs <- data.frame(
  group = "Legumes",
  n = N,
  lower_bound = lower_bound,
  mid_point = best_est,
  upper_bound = upper_bound,
  dataset = "Predictions",
  year = "2024"
)

# monocots
mons_preds <- preds_raw %>%
  left_join(APGIV_families, by = c(family = 'APG.IV.Family')) %>%
  filter(HigherGroups == 'Monocots', complete.cases(combined)) %>%
  select(family, threatened, .lower, .upper, confidence, combined)

mons_summary <- mons_preds %>% group_by(combined, confidence) %>% count()

N <- sum(mons_summary$n)

# Extract counts
threat_conf  <- mons_summary$n[mons_summary$combined == "threatened" &
                                 mons_summary$confidence == "confident"]

threat_low<- mons_summary$n[mons_summary$combined == "threatened" &
                              mons_summary$confidence == "low_confidence"]

not_conf     <- mons_summary$n[mons_summary$combined == "not threatened" &
                                 mons_summary$confidence == "confident"]

# Calculations
lower_bound  <- threat_conf / N
upper_bound  <- (threat_conf + threat_low) / N
best_est     <- threat_conf / (threat_conf + not_conf)

# Combine results
RL_2024_preds_perc_threat_mons <- data.frame(
  group = "Monocots",
  n = N,
  lower_bound = lower_bound,
  mid_point = best_est,
  upper_bound = upper_bound,
  dataset = "Predictions",
  year = "2024"
)

RL_2024_preds_perc_threat <- bind_rows(RL_2024_preds_perc_threat_mons,
                                       RL_2024_preds_perc_threat_legs)

# Sup Tables - categories  ----
RL_categories_per_group <- bind_rows(
  SRLI_2020_categories, # SRLI 2020
  SRLI_2010_categories, # SRLI original
  RL_2010_categories, # RL 2010
  RL_2015_categories, # RL 2015
  RL_2020_categories, # RL 2020
  RL_2024_categories, # RL 2024
  RL_2024_preds_categories # RL 2024 + preds
)

# add percentage
RL_categories_per_group <- RL_categories_per_group %>%
  group_by(dataset, year, group) %>%
  mutate(
    total = sum(count, na.rm = TRUE),
    percent = count / total * 100
  ) %>%
  ungroup()  %>%
  mutate(percent = round(percent, 2))

# write this table to results and merge as a single excel file in sup mat
write.csv(RL_categories_per_group,
          "02_figures/Table_S2_RL_categories_per_group.csv",
          row.names = FALSE)

# stat tests ----
# compare DDs between groups
test_data <- RL_categories_per_group %>%
  group_by(group) %>%
  filter(year == 2000) %>%
  filter(category != "NE") %>%
  summarise(
    DD = sum(count[category == "DD"], na.rm = TRUE),
    total = sum(count, na.rm = TRUE),
    non_DD = total - DD,
    .groups = "drop"
  ) %>%
  mutate(DD_perc = round((DD / total) * 100, 1))

mat <- as.matrix(test_data[, c("DD", "non_DD")])
rownames(mat) <- test_data$group

chisq.test(mat)

legs_perc_DD <-

# Sup Tables - percentage threatened  ----
RL_perc_threat <- bind_rows(
  SRLI_2020_perc_threat, # SRLI 2020
  #SRLI_2020_perc_threat_augmented, # SRLI 2020 + preds
  SRLI_2010_perc_threat, # SRLI original
  RL_2010_perc_threat, # RL 2010
  RL_2015_perc_threat, # RL 2010
  RL_2020_perc_threat, # RL 2020
  RL_2024_perc_threat, # RL 2024
  RL_2024_preds_perc_threat
  #RL_2024_DD_perc_threat
  #Preds_2024_perc_threat # preds 2024
)

# write this table to results and merge as a single excel file in sup mat
write.csv(RL_perc_threat,
          "02_figures/Table_S3_percentage_threatened.csv",
          row.names = FALSE)

# Figure threat percentages ----
# add n for each dataset

# set this to get the order correct
RL_perc_threat <- RL_perc_threat %>%
  mutate(dataset = factor(
    dataset,
    levels = c("SRLI 2010", "SRLI 2020", "Red List","Predictions") #"Red List 2024 +NE"
  )) %>%
  mutate(n_scaled = rescale(n, to = c(0, 1)))

pd <- position_dodge(width = 1)

ggplot(RL_perc_threat,
       aes(x = as.numeric(year),
           y = mid_point,
           colour = group,
           group = interaction(group, dataset))) +
  facet_wrap(~ dataset, nrow = 1, labeller = label_wrap_gen(width = 12)) +

  geom_line(position = pd) +
  # wide uncertainty (drawn first → behind)
  # geom_linerange(
  #   aes(ymin = NE_lower, ymax = NE_upper),
  #   position = pd,
  #   colour = "grey70",
  #   linewidth = 0.6
  # ) +
  geom_pointrange(
    aes(ymin = lower_bound, ymax = upper_bound, colour = group),
    position = pd
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format()
  ) +
  scale_color_manual(values = c(
    "Monocots" = "#0072B2",
    "Legumes"  = "#E69F00"
  ), name = NULL) +
  labs(x = "Year", y = "Proportion threatened") +
  theme_base(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    strip.background = element_rect(fill = "white", colour = NA),
    strip.text = element_text(face = "bold")
  ) +
  scale_x_continuous(n.breaks = 3)

ggsave("02_figures/Figure_1_proportion_threatened.png",
        width = 8,
        height = 4)




