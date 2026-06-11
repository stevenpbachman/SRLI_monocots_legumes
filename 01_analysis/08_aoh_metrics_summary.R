

# AoH metrics for genuine changes ----

# libraries -
library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(ggbeeswarm)
library(ggdist)
library(ggExtra)
library(hexbin)
library(ggthemes)
library(stringr)
library(purrr)
library(rCAT)

# Summary of AOH changes to identify genuine changes

# organise data
# start from scratch so that this is standalone
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir, "00_raw_data/SRLI_2026.csv")) %>%
  dplyr::filter(SRLI_2024_analysis == "include")

# load functions below #results summary etc.
# single to test it is working
powo_ID = "469774-1"
sp_res <- read_csv(paste0(work_dir,"04_outputs/curated_workflow/aoh_stats/results_", powo_ID ,".csv"))
sp_res <- read_csv(paste0(work_dir,"/04_outputs/gbif_workflow/aoh_stats/results_", powo_ID ,".csv"))

# batch generate results
# Read in the AOH results
# GBIF first ----
files <- list.files(path = paste0(work_dir, "04_outputs/gbif_workflow/aoh_stats/"),
                    pattern="csv",
                    all.files=FALSE,
                    full.names=TRUE,
                    recursive=TRUE)
thefile = str_match(files, "results_(.*?).csv")
results = thefile[,2]
length(results)

to_do <- results

i = 1
sp_res <- read_csv(paste0(work_dir,"04_outputs/gbif_workflow/aoh_stats/results_", to_do[i] ,".csv"))
res = results_summary(sp_res)

for (i in 2:length(to_do)){

  sp_res <- read_csv(paste0(work_dir,"04_outputs/gbif_workflow/aoh_stats/results_", to_do[i] ,".csv"),
                     col_types = cols()) # to suppress warnings
  print(sp_res$powoID[i])
  res_l = try(results_summary(sp_res)) #silent=TRUE
  res = bind_rows(res, res_l)

}

# save
gbif_res_all <- res
#write_csv(gbif_res_all, paste0("04_outputs/all_gbif_possible_genuine_changes.csv"))

# work out only the genuine possible changes i.e. LC, NT cancelled out by VU

# filter to get the possible changes
res_comb_gbif <- res %>%
  dplyr::filter(gen_aoh_2000_2010_10yr == TRUE |
                  gen_aoh_2010_2020_10yr == TRUE |
                  gen_aoh_2000_2010_20yr == TRUE |
                  gen_aoh_2010_2020_20yr == TRUE |
                  gen_aoh_2000_2010_30yr == TRUE |
                  gen_aoh_2010_2020_30yr == TRUE |
                  gen_aoo_max_2000_2010== TRUE |
                  gen_aoo_max_2010_2020 == TRUE |
                  gen_eoo_max_2000_2010 == TRUE |
                  gen_eoo_max_2010_2020 == TRUE
                ) %>%
  mutate(workflow = "gbif")


# batch generate results
# Read in the AOH results
# Curated second ----
files <- list.files(path = paste0(work_dir, "04_outputs/curated_workflow/aoh_stats/"),
                    pattern="csv",
                    all.files=FALSE,
                    full.names=TRUE,
                    recursive=TRUE)
thefile = str_match(files, "results_(.*?).csv")
results = thefile[,2]
length(results)

to_do <- results

i = 1
sp_res <- read_csv(paste0(work_dir,"04_outputs/curated_workflow/aoh_stats/results_", to_do[i] ,".csv"))
res = results_summary(sp_res)

for (i in 2:length(to_do)){

  sp_res <- read_csv(paste0(work_dir,"04_outputs/curated_workflow/aoh_stats/results_", to_do[i] ,".csv"),
                     col_types = cols()) # to suppress warnings
  print(sp_res$powoID[i])
  res_l = try(results_summary(sp_res)) #silent=TRUE
  res = bind_rows(res, res_l)

}

# save
curated_res_all <- res

# work out only the genuine possible changes i.e. LC, NT cancelled out by VU

# filter to get the possible changes
res_comb_curated <- res %>%
  dplyr::filter(gen_aoh_2000_2010_10yr == TRUE |
                  gen_aoh_2010_2020_10yr == TRUE |
                  gen_aoh_2000_2010_20yr == TRUE |
                  gen_aoh_2010_2020_20yr == TRUE |
                  gen_aoh_2000_2010_30yr == TRUE |
                  gen_aoh_2010_2020_30yr == TRUE |
                  gen_aoo_max_2000_2010== TRUE |
                  gen_aoo_max_2010_2020 == TRUE |
                  gen_eoo_max_2000_2010 == TRUE |
                  gen_eoo_max_2010_2020 == TRUE
                ) %>%
  mutate(workflow = "curated")


# summarise/combine  results
res_comb <- bind_rows(res_comb_gbif, res_comb_curated)

res_comb_true <- res_comb %>% dplyr::select("gen_eoo_max_2000_2010",
                                            "gen_eoo_max_2010_2020",
                                            "gen_aoo_max_2000_2010",
                                            "gen_aoo_max_2010_2020",
                                            "gen_aoh_2000_2010_30yr",
                                            "gen_aoh_2010_2020_30yr",
                                            "gen_aoh_2000_2010_20yr",
                                            "gen_aoh_2010_2020_20yr",
                                            "gen_aoh_2000_2010_10yr",
                                            "gen_aoh_2010_2020_10yr"
)



# save as a table?
# Summarise ----
sum_res <- data.frame(
  column = names(res_comb_true),
  n_true = colSums(res_comb_true, na.rm = TRUE)
)

# how many species?
res_comb_sp <- res_comb %>% group_by() %>% distinct(powo_ID) %>% count()

# # join with SRLI to get 2020_ratings
# SRLI_2024_thin <- SRLI_2024 %>% select(WCVP_current_id, RL_2000, RL_2010, RL_2020, WCVP_accepted_name)
# res_comb_id <- res_comb %>% left_join(SRLI_2024_thin, by = c("powo_ID" = "WCVP_current_id"))
# test <- res_comb_id %>% group_by(powo_ID) %>% count()
# write_csv(res_comb_id, paste0("G:/Users/sb42kg/SRLI_2024/04_outputs/possible_genuine_changes_curated.csv"))
# write_csv(res_comb_id, paste0("G:/Users/sb42kg/SRLI_2024/04_outputs/possible_genuine_changes_gbif.csv"))
# #write_csv(res_comb_id, paste0("G:/Users/sb42kg/SRLI_2024/04_outputs/possible_genuine_changes.csv"))






#### functions ####
results_summary <- function(sp_res){

  powo_ID <- sp_res$powoID[1]

  # prepare the SHI fields
  SHI_aoh_cells_2020 = sp_res$cell_count[29]
  aoh_prop_full = sp_res$prop_change[29]

  # summarise the proportion change statistics
  # 10 year decline at 3 time points
  # old way
  #aoh_prop_2000_10yr = (sp_res$prop_change[9] -  sp_res$prop_change[1]) * 100
  #aoh_prop_2010_10yr = (sp_res$prop_change[19] -  sp_res$prop_change[10]) * 100
  #aoh_prop_2020_10yr = (sp_res$prop_change[29] -  sp_res$prop_change[20])* 100

  # new way
  aoh_prop_2000_10yr = ((sp_res$cell_count[9] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100
  aoh_prop_2010_10yr = ((sp_res$cell_count[19] -  sp_res$cell_count[10]) / sp_res$cell_count[10])*100
  aoh_prop_2020_10yr = ((sp_res$cell_count[29] -  sp_res$cell_count[20]) / sp_res$cell_count[20])*100

  aoh_2000_10yr_pop_trend = pop_trend(aoh_prop_2000_10yr)
  aoh_2010_10yr_pop_trend = pop_trend(aoh_prop_2010_10yr)
  aoh_2020_10yr_pop_trend = pop_trend(aoh_prop_2020_10yr)

  # get max category
  if (aoh_2000_10yr_pop_trend == "Decreasing"){ aoh_trend_cat_2000_10yr = rCAT::ratingPop(abs(aoh_prop_2000_10yr))
  } else {aoh_trend_cat_2000_10yr = "LC"}

  if (aoh_2010_10yr_pop_trend == "Decreasing"){ aoh_trend_cat_2010_10yr = rCAT::ratingPop(abs(aoh_prop_2010_10yr))
  } else {aoh_trend_cat_2010_10yr = "LC"}

  if (aoh_2020_10yr_pop_trend == "Decreasing"){aoh_trend_cat_2020_10yr = rCAT::ratingPop(abs(aoh_prop_2020_10yr))
  } else {aoh_trend_cat_2020_10yr = "LC"}

  gen_aoh_2000_2010_10yr <- aoh_trend_cat_2000_10yr != aoh_trend_cat_2010_10yr
  gen_aoh_2010_2020_10yr <- aoh_trend_cat_2010_10yr != aoh_trend_cat_2020_10yr

  # 20 year decline at 3 time points
  aoh_prop_2000_20yr = ((sp_res$cell_count[9] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100
  aoh_prop_2010_20yr = ((sp_res$cell_count[19] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100
  aoh_prop_2020_20yr = ((sp_res$cell_count[29] -  sp_res$cell_count[10])/ sp_res$cell_count[10])*100

  aoh_2000_20yr_pop_trend = pop_trend(aoh_prop_2000_20yr)
  aoh_2010_20yr_pop_trend = pop_trend(aoh_prop_2010_20yr)
  aoh_2020_20yr_pop_trend = pop_trend(aoh_prop_2020_20yr)

  # get max category
  if (aoh_2000_20yr_pop_trend == "Decreasing"){ aoh_trend_cat_2000_20yr = rCAT::ratingPop(abs(aoh_prop_2000_20yr))
  } else {aoh_trend_cat_2000_20yr = "LC"}

  if (aoh_2010_20yr_pop_trend == "Decreasing"){ aoh_trend_cat_2010_20yr = rCAT::ratingPop(abs(aoh_prop_2010_20yr))
  } else {aoh_trend_cat_2010_20yr = "LC"}

  if (aoh_2020_20yr_pop_trend == "Decreasing"){aoh_trend_cat_2020_20yr = rCAT::ratingPop(abs(aoh_prop_2020_20yr))
  } else {aoh_trend_cat_2020_20yr = "LC"}

  gen_aoh_2000_2010_20yr <- aoh_trend_cat_2000_20yr != aoh_trend_cat_2010_20yr
  gen_aoh_2010_2020_20yr <- aoh_trend_cat_2010_20yr != aoh_trend_cat_2020_20yr

  # 30 year decline at 3 time points
  aoh_prop_2000_30yr = ((sp_res$cell_count[9] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100
  aoh_prop_2010_30yr = ((sp_res$cell_count[19] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100
  aoh_prop_2020_30yr = ((sp_res$cell_count[29] -  sp_res$cell_count[1]) / sp_res$cell_count[1])*100

  aoh_2000_30yr_pop_trend = pop_trend(aoh_prop_2000_30yr)
  aoh_2010_30yr_pop_trend = pop_trend(aoh_prop_2010_30yr)
  aoh_2020_30yr_pop_trend = pop_trend(aoh_prop_2020_30yr)

  # get max category
  if (aoh_2000_30yr_pop_trend == "Decreasing"){ aoh_trend_cat_2000_30yr = rCAT::ratingPop(abs(aoh_prop_2000_30yr))
  } else {aoh_trend_cat_2000_30yr = "LC"}

  if (aoh_2010_30yr_pop_trend == "Decreasing"){ aoh_trend_cat_2010_30yr = rCAT::ratingPop(abs(aoh_prop_2010_30yr))
  } else {aoh_trend_cat_2010_30yr = "LC"}

  if (aoh_2020_30yr_pop_trend == "Decreasing"){aoh_trend_cat_2020_30yr = rCAT::ratingPop(abs(aoh_prop_2020_30yr))
  } else {aoh_trend_cat_2020_30yr = "LC"}

  gen_aoh_2000_2010_30yr <- aoh_trend_cat_2000_30yr != aoh_trend_cat_2010_30yr
  gen_aoh_2010_2020_30yr <- aoh_trend_cat_2010_30yr != aoh_trend_cat_2020_30yr


  # Summarise EOO results
  #eoo_mcp_occs_2000 <- sp_res$eoo_aoh[9]
  #eoo_mcp_occs_2010 <- sp_res$eoo_aoh[19]
  #eoo_mcp_occs_2020 <- sp_res$eoo_aoh[29]

  eoo_aoh_2000 <- sp_res$eoo_aoh[9]
  eoo_aoh_2010 <- sp_res$eoo_aoh[19]
  eoo_aoh_2020 <- sp_res$eoo_aoh[29]

  # 1. check for genuine change using the max method - evidentiary
  eoo_max_cat_2000 <- rCAT::ratingEoo(eoo_aoh_2000)
  eoo_max_cat_2010 <- rCAT::ratingEoo(eoo_aoh_2010)
  eoo_max_cat_2020 <- rCAT::ratingEoo(eoo_aoh_2020)

  gen_eoo_max_2000_2010 <- eoo_max_cat_2000 != eoo_max_cat_2010
  gen_eoo_max_2010_2020 <- eoo_max_cat_2010 != eoo_max_cat_2020

  # 2. check for genuine change using the majority method - precautionary
  # eoo_range_2000 <- round(min(eoo_mcp_buff_2000, eoo_mcp_occs_2000)):round(max(eoo_mcp_occs_2000, eoo_mcp_buff_2000))
  # eoo_maj_cat_2000 <- majority_cat_eoo(eoo_range_2000)
  #
  # eoo_range_2010 <- round(min(eoo_mcp_buff_2010, eoo_mcp_occs_2010)):round(max(eoo_mcp_occs_2010, eoo_mcp_buff_2010))
  # eoo_maj_cat_2010 <- majority_cat_eoo(eoo_range_2010)
  #
  # eoo_range_2020 <- round(min(eoo_mcp_buff_2020, eoo_mcp_occs_2020)):round(max(eoo_mcp_occs_2020, eoo_mcp_buff_2020))
  # eoo_maj_cat_2020 <- majority_cat_eoo(eoo_range_2020)
  #
  # gen_eoo_maj_2000_2010 <- eoo_maj_cat_2000 != eoo_maj_cat_2010
  # gen_eoo_maj_2010_2020 <- eoo_maj_cat_2010 != eoo_maj_cat_2020

  # 3. check for genuine change using the average(mean) method
  # eoo_mean_2000 <- mean(c(eoo_mcp_occs_2000,eoo_mcp_buff_2000))
  # eoo_mean_cat_2000 <- rCAT::ratingEoo(eoo_mean_2000)
  #
  # eoo_mean_2010 <- mean(c(eoo_mcp_occs_2010,eoo_mcp_buff_2010))
  # eoo_mean_cat_2010 <- rCAT::ratingEoo(eoo_mean_2010)
  #
  # eoo_mean_2020 <- mean(c(eoo_mcp_occs_2020,eoo_mcp_buff_2020))
  # eoo_mean_cat_2020 <- rCAT::ratingEoo(eoo_mean_2010)
  #
  # gen_eoo_mean_2000_2010 <- eoo_mean_cat_2000 != eoo_mean_cat_2010
  # gen_eoo_mean_2010_2020 <- eoo_mean_cat_2010 != eoo_mean_cat_2020

  # summarise AOO results
  # aoo_occs_2000 <- sp_res$aoo_occs[9]
  # aoo_occs_2010 <- sp_res$aoo_occs[19]
  # aoo_occs_2020 <- sp_res$aoo_occs[29]

  aoo_aoh_2000 <- sp_res$aoo_aoh[9]
  aoo_aoh_2010 <- sp_res$aoo_aoh[19]
  aoo_aoh_2020 <- sp_res$aoo_aoh[29]

  # 1. check for genuine change using the max method - evidentiary
  # aoo_max_cat_2000 <- rCAT::ratingAoo(max(aoo_occs_2000,aoo_aoh_2000))
  # aoo_max_cat_2010 <- rCAT::ratingAoo(max(aoo_occs_2010,aoo_aoh_2010))
  # aoo_max_cat_2020 <- rCAT::ratingAoo(max(aoo_occs_2020,aoo_aoh_2020))

  aoo_max_cat_2000 <- rCAT::ratingAoo(aoo_aoh_2000)
  aoo_max_cat_2010 <- rCAT::ratingAoo(aoo_aoh_2010)
  aoo_max_cat_2020 <- rCAT::ratingAoo(aoo_aoh_2020)

  gen_aoo_max_2000_2010 <- aoo_max_cat_2000 != aoo_max_cat_2010
  gen_aoo_max_2010_2020 <- aoo_max_cat_2010 != aoo_max_cat_2020


  # 2. check for genuine change using the majority method - precautionary
  # aoo_range_2000 <- round(min(aoo_aoh_2000, aoo_occs_2000)):round(max(aoo_occs_2000, aoo_aoh_2000))
  # aoo_maj_cat_2000 <- majority_cat_aoo(aoo_range_2000)
  #
  # aoo_range_2010 <- round(min(aoo_aoh_2010, aoo_occs_2010)):round(max(aoo_occs_2010, aoo_aoh_2010))
  # aoo_maj_cat_2010 <- majority_cat_aoo(aoo_range_2010)
  #
  # aoo_range_2020 <- round(min(aoo_aoh_2020, aoo_occs_2020)):round(max(aoo_occs_2020, aoo_aoh_2020))
  # aoo_maj_cat_2020 <- majority_cat_aoo(aoo_range_2020)
  #
  # gen_aoo_maj_2000_2010 <- aoo_maj_cat_2000 != aoo_maj_cat_2010
  # gen_aoo_maj_2010_2020 <- aoo_maj_cat_2010 != aoo_maj_cat_2020

  # 3. check for genuine change using the average(mean) method
  # aoo_mean_2000 <- mean(c(aoo_occs_2000,aoo_aoh_2000))
  # aoo_mean_cat_2000 <- rCAT::ratingAoo(aoo_mean_2000)
  #
  # aoo_mean_2010 <- mean(c(aoo_occs_2010,aoo_aoh_2010))
  # aoo_mean_cat_2010 <- rCAT::ratingAoo(aoo_mean_2010)
  #
  # aoo_mean_2020 <- mean(c(aoo_occs_2020,aoo_aoh_2020))
  # aoo_mean_cat_2020 <- rCAT::ratingAoo(aoo_mean_2020)
  #
  # gen_aoo_mean_2000_2010 <- aoo_mean_cat_2000 != aoo_mean_cat_2010
  # gen_aoo_mean_2010_2020 <- aoo_mean_cat_2010 != aoo_mean_cat_2020

  results <- data.frame(
    "powo_ID" = powo_ID,
    "aoh_prop_2000_10yr" = aoh_prop_2000_10yr,
    "aoh_prop_2010_10yr" = aoh_prop_2010_10yr,
    "aoh_prop_2020_10yr" = aoh_prop_2020_10yr,
    "aoh_2000_pop_trend_10yr" = aoh_2000_10yr_pop_trend,
    "aoh_2010_pop_trend_10yr" = aoh_2010_10yr_pop_trend,
    "aoh_2020_pop_trend_10yr" = aoh_2020_10yr_pop_trend,
    "aoh_trend_cat_2000_10yr" = aoh_trend_cat_2000_10yr,
    "aoh_trend_cat_2010_10yr" = aoh_trend_cat_2010_10yr,
    "aoh_trend_cat_2020_10yr" = aoh_trend_cat_2020_10yr,
    "gen_aoh_2000_2010_10yr" = gen_aoh_2000_2010_10yr,
    "gen_aoh_2010_2020_10yr" = gen_aoh_2010_2020_10yr,

    "aoh_prop_2000_20yr" = aoh_prop_2000_20yr,
    "aoh_prop_2010_20yr" = aoh_prop_2010_20yr,
    "aoh_prop_2020_20yr" = aoh_prop_2020_20yr,
    "aoh_2000_pop_trend_20yr" = aoh_2000_20yr_pop_trend,
    "aoh_2010_pop_trend_20yr" = aoh_2010_20yr_pop_trend,
    "aoh_2020_pop_trend_20yr" = aoh_2020_20yr_pop_trend,
    "aoh_trend_cat_2000_20yr" = aoh_trend_cat_2000_20yr,
    "aoh_trend_cat_2010_20yr" = aoh_trend_cat_2010_20yr,
    "aoh_trend_cat_2020_20yr" = aoh_trend_cat_2020_20yr,
    "gen_aoh_2000_2010_20yr" = gen_aoh_2000_2010_20yr,
    "gen_aoh_2010_2020_20yr" = gen_aoh_2010_2020_20yr,

    "aoh_prop_2000_30yr" = aoh_prop_2000_30yr,
    "aoh_prop_2010_30yr" = aoh_prop_2010_30yr,
    "aoh_prop_2020_30yr" = aoh_prop_2020_30yr,
    "aoh_2000_pop_trend_30yr" = aoh_2000_30yr_pop_trend,
    "aoh_2010_pop_trend_30yr" = aoh_2010_30yr_pop_trend,
    "aoh_2020_pop_trend_30yr" = aoh_2020_30yr_pop_trend,
    "aoh_trend_cat_2000_30yr" = aoh_trend_cat_2000_30yr,
    "aoh_trend_cat_2010_30yr" = aoh_trend_cat_2010_30yr,
    "aoh_trend_cat_2020_30yr" = aoh_trend_cat_2020_30yr,
    "gen_aoh_2000_2010_30yr" = gen_aoh_2000_2010_30yr,
    "gen_aoh_2010_2020_30yr" = gen_aoh_2010_2020_30yr,

    "aoo_aoh_2000"=aoo_aoh_2000,
    "aoo_aoh_2010"=aoo_aoh_2010,
    "aoo_aoh_2020"=aoo_aoh_2020,
    # "aoo_mean_cat_2000" = aoo_mean_cat_2000,
    # "aoo_mean_cat_2010" = aoo_mean_cat_2010,
    # "aoo_mean_cat_2020" = aoo_mean_cat_2020,
    # "gen_aoo_mean_2000_2010" = gen_aoo_mean_2000_2010,
    # "gen_aoo_mean_2010_2020" = gen_aoo_mean_2010_2020,
    "aoo_max_cat_2000"=aoo_max_cat_2000,
    "aoo_max_cat_2010"=aoo_max_cat_2010,
    "aoo_max_cat_2020"=aoo_max_cat_2020,
    "gen_aoo_max_2000_2010"=gen_aoo_max_2000_2010,
    "gen_aoo_max_2010_2020"=gen_aoo_max_2010_2020,
    # "aoo_maj_cat_2000"=aoo_maj_cat_2000,
    # "aoo_maj_cat_2010"=aoo_maj_cat_2010,
    # "aoo_maj_cat_2020"=aoo_maj_cat_2020,
    # "gen_aoo_maj_2000_2010"=gen_aoo_maj_2000_2010,
    # "gen_aoo_maj_2010_2020"=gen_aoo_maj_2010_2020,
    # "eoo_mean_cat_2000" = eoo_mean_cat_2000,
    # "eoo_mean_cat_2010" = eoo_mean_cat_2010,
    # "eoo_mean_cat_2020" = eoo_mean_cat_2020,
    # "gen_eoo_mean_2000_2010" = gen_eoo_mean_2000_2010,
    # "gen_eoo_mean_2010_2020" = gen_eoo_mean_2010_2020,
    "eoo_max_cat_2000" = eoo_max_cat_2000,
    "eoo_max_cat_2010" = eoo_max_cat_2010,
    "eoo_max_cat_2020" = eoo_max_cat_2020,
    "gen_eoo_max_2000_2010"=gen_eoo_max_2000_2010,
    "gen_eoo_max_2010_2020"=gen_eoo_max_2010_2020,
    "eoo_aoh_2000" = eoo_aoh_2000,
    "eoo_aoh_2010" = eoo_aoh_2010,
    "eoo_aoh_2020" = eoo_aoh_2020
    # "eoo_maj_cat_2000"=eoo_maj_cat_2000,
    # "eoo_maj_cat_2010"=eoo_maj_cat_2010,
    # "eoo_maj_cat_2020"=eoo_maj_cat_2020,
    # "gen_eoo_maj_2000_2010"=gen_eoo_maj_2000_2010,
    # "gen_eoo_maj_2010_2020"=gen_eoo_maj_2010_2020
    )

  return(results)


}

#range_value = aoo_range_2000

majority_cat_eoo <- function(range_value){

  range_value <- as.numeric(range_value)
  cats <- map(range_value, rCAT::ratingEoo)
 # cats <- lapply(range_value,rCAT::ratingEoo)
  cats_df <- as.data.frame(do.call(rbind, cats)) %>%
    group_by(V1) %>%
    summarise(Percentage=n()/nrow(.))

  maj_cat <- cats_df %>%
    slice(which.max(Percentage)) %>%
    select(V1) %>%
    as.character()

  return(maj_cat)

}
majority_cat_aoo <- function(range_value){

  range_value <- as.numeric(range_value)
  cats <- map(range_value, rCAT::ratingAoo)
  #cats <- lapply(range_value,rCAT::ratingAoo)
    cats_df <- as.data.frame(do.call(rbind, cats)) %>%
    group_by(V1) %>%
    summarise(Percentage=n()/nrow(.))

  maj_cat <- cats_df %>%
    slice(which.max(Percentage)) %>%
    select(V1) %>%
    as.character()

  return(maj_cat)

}
pop_trend <-function(pop_data){

  # plus or minus 1% = stable
  # >1 = increasing
  # <1 = decreasing

    if (pop_data > 1){
      cat <- "Increasing"
    } else if (pop_data < -1){
      cat <- "Decreasing"
    } else
      cat <- "Stable"

}



## explore possible chart to show majority change and maximum change
#ggplot(cats_df, aes(fill=condition, y=value, x=specie)) +
#  geom_bar(position="dodge", stat="identity")

################## EXTRAS ###############

#EOO_range = filter on eoo_mcp_occs_cat & eoo_mcp_occs_cat. Use unique if they are the same category
aoo_range = sort(c(sp_res$aoo_occs[29], sp_res$aoo_buff_aoh[29]))
aoo_range = paste0(aoo_range[1],"-", aoo_range[2])
aoo_range_cat = aoo_range_cat(sp_res)

eoo_range = sort(c(sp_res$eoo_mcp_occs[29], sp_res$eoo_mcp_buff[29]))
eoo_range = paste0(eoo_range[1],"-", eoo_range[2])
eoo_range_cat = eoo_range_cat(sp_res)


aoo_range_cat = function(sp_res){

  aoo_range<- sp_res %>%
  dplyr::filter(year == "2020") %>%
  select(aoo_occs_cat,aoo_buff_aoh_cat)

occs <- case_when(
  aoo_range$aoo_occs_cat == "LC" ~ "5",
  aoo_range$aoo_occs_cat == "NT" ~ "4",
  aoo_range$aoo_occs_cat == "VU" ~ "3",
  aoo_range$aoo_occs_cat == "EN" ~ "2",
  aoo_range$aoo_occs_cat == "CR" ~ "1"
)

buff <- case_when(
  aoo_range$aoo_buff_aoh_cat == "LC" ~ "5",
  aoo_range$aoo_buff_aoh_cat == "NT" ~ "4",
  aoo_range$aoo_buff_aoh_cat == "VU" ~ "3",
  aoo_range$aoo_buff_aoh_cat == "EN" ~ "2",
  aoo_range$aoo_buff_aoh_cat == "CR" ~ "1"
)

t = sort(c(buff, occs))

r <- case_when(
  t == 5 ~ "LC",
  t == 4 ~ "NT",
  t == 3 ~ "VU",
  t == 2 ~ "EN",
  t == 1 ~ "CR"
)}
eoo_range_cat = function(sp_res){

  eoo_range<- sp_res %>%
    dplyr::filter(year == "2020") %>%
    select(eoo_mcp_occs_cat,eoo_mcp_buff_cat)

  occs <- case_when(
    eoo_range$eoo_mcp_occs_cat  == "LC" ~ "5",
    eoo_range$eoo_mcp_occs_cat == "NT" ~ "4",
    eoo_range$eoo_mcp_occs_cat == "VU" ~ "3",
    eoo_range$eoo_mcp_occs_cat == "EN" ~ "2",
    eoo_range$eoo_mcp_occs_cat == "CR" ~ "1"
  )

  buff <- case_when(
    eoo_range$eoo_mcp_buff_cat == "LC" ~ "5",
    eoo_range$eoo_mcp_buff_cat == "NT" ~ "4",
    eoo_range$eoo_mcp_buff_cat == "VU" ~ "3",
    eoo_range$eoo_mcp_buff_cat == "EN" ~ "2",
    eoo_range$eoo_mcp_buff_cat == "CR" ~ "1"
  )

  t = sort(c(buff, occs))

  r <- case_when(
    t == 5 ~ "LC",
    t == 4 ~ "NT",
    t == 3 ~ "VU",
    t == 2 ~ "EN",
    t == 1 ~ "CR"
  )}




