
# Make figures and save outputs in:
# SRLI_project\SRLI_2_paper_repo\02_figures

# Simple script to generate all figures
# Use this to control the order

library(writexl)

# source the figures scripts
source("01_analysis/01_asynchronous_assessments.R") # Figure S1
source("01_analysis/02_data_coverage.R") # Table S1, Figure S2
source("01_analysis/03_threat_proportion.R") # Figure 1, Tables S2,S3
source("01_analysis/04_category_changes.R") # Table 1, Figure 2 (A,B,C,D), Tables S4,S5
source("01_analysis/05_RLI_trend_lines.R") #Figure 3, Tables S6
source("01_analysis/06_compare_curated_gbif_aoh.R") #Figure S3
source("01_analysis/07_aoh_validation_results.R") #Figure 4, Table_S7, Figure S4


# generate the supplementary tables - excel
data_coverage_summary <- read.csv("02_figures/Table_S1_data_coverage_summary.csv")
RL_categories_per_group <- read.csv("02_figures/Table_S2_RL_categories_per_group.csv")
RL_perc_threat <- read.csv("02_figures/Table_S3_percentage_threatened.csv")
change_summary <- read.csv("02_figures/change_summary.csv")
genuine_threat_reasons <- read.csv("02_figures/Table_S4_genuine_threat_reasons.csv")
non_genuine_raw_counts <- read.csv("02_figures/Table_S5_non_genuine_raw_counts.csv")
RLI_values <- read.csv("02_figures/Table_S6_rli_values.csv")
Validation_proportions <- read.csv("02_figures/Table_S7_validation_proportions.csv")

# write as an excel file with tabs for each table
write_xlsx(
  list("S1_data_coverage_summary" = data_coverage_summary,
  "S2_RL_categories_per_group" = RL_categories_per_group,
  "S3_RL_perc_threat" = RL_perc_threat,
  "S4_genuine_threat_reasons" = genuine_threat_reasons,
  "S5_non_genuine_raw_counts.csv" = non_genuine_raw_counts,
  "S6_RLI_values.csv" = RLI_values,
  "S7_Validation_proportions.csv" = Validation_proportions
  ),
  "02_figures/Supplementary_tables.xlsx")
