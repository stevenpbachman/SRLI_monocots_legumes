
# Figure S3 ----

# Comparison of range size of species by data source (GBIF vs curated), as
# measured by AoH cell count in 2020.

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

# organise data
# start from scratch so that this is standalone
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir, "00_raw_data/SRLI_2026.csv")) %>%
  dplyr::filter(SRLI_2024_analysis == "include")

# use the data from part 1. to feed part 2

# 1. Basic stats on number of AOH maps ----
# GBIF workflow ----
# total number of AOH maps
aoh_files_gbif <- list.files(path=paste0(work_dir,"04_outputs/gbif_workflow/aoh_stats/"),
                             pattern="csv",
                             all.files=FALSE,
                             full.names=TRUE,
                             recursive=TRUE)

# read in the AOH summary files
gbif_aoh_read_in = purrr::map_df(aoh_files_gbif, read_csv, show_col_types = FALSE)

# save combined file for later
aoh_gbif <- data.frame(powoID = unique(gbif_aoh_read_in$powoID), aoh_GBIF_stats = "Yes")

#test_dups <- aoh_gbif %>% group_by(powoID) %>% count()
print(paste0("AOH result 1: ", nrow(aoh_gbif), " GBIF AOH maps"))

# curated workflow ----
# total number of AOH maps
aoh_files_curated <- list.files(path=paste0(work_dir,"04_outputs/curated_workflow/aoh_stats/"),
                                pattern="csv",
                                all.files=FALSE,
                                full.names=TRUE,
                                recursive=TRUE)

# read in the AOH summary files
curated_aoh_read_in = purrr::map_df(aoh_files_curated, read_csv, show_col_types = FALSE) %>%
  dplyr::mutate(source = "curated")

# save combined file for later
aoh_curated <- data.frame(powoID = unique(curated_aoh_read_in$powoID), aoh_curated_stats = "Yes")

test_dups <- aoh_curated %>% group_by(powoID) %>% count()
print(paste0("AOH result 1: ", nrow(aoh_curated), " curated AOH maps"))

# join them together and summarise
aoh_comb <- aoh_gbif %>% left_join(aoh_curated, by = 'powoID')

# go long to get total number of aoh maps
aoh_comb_long <- bind_rows(aoh_gbif, aoh_curated)
print(paste0("AOH result 2: ", nrow(aoh_comb_long), " AOH maps"))

# join to get counts per groups
comb_long <- aoh_comb_long %>% left_join(SRLI_2026,
                                            by = c("powoID" = "WCVP_current_id")) %>%
  select(group, powoID)

count_groups <- comb_long %>%
  group_by(group) %>% count()


# how many species in both
print(paste0("AOH result 3: ",
             sum(aoh_comb$aoh_GBIF_stats == "Yes" &
                   aoh_comb$aoh_curated_stats == "Yes",
                 na.rm = TRUE),
             " both GBIF and curated AOH maps"))


# 2. GBIF vs curated validation ----

curated_2020 <- curated_aoh_read_in %>% filter(year == "year_2020") %>%
  select(cell_count, powoID) %>%
  mutate(workflow = "curated")

GBIF_2020 <- gbif_aoh_read_in %>% filter(year == "year_2020") %>%
  select(cell_count, powoID) %>%
  mutate(workflow = "gbif")

comb_2020 <- bind_rows(curated_2020, GBIF_2020)

# join to get counts per groups
comb_2020_groups <- comb_2020 %>% left_join(SRLI_2026,
                                            by = c("powoID" = "WCVP_current_id")) %>%
select(group, cell_count, powoID, workflow)

count_groups <- comb_2020_groups %>%
  group_by(group) %>% count()

# Keep only powoIDs that have both 'curated' and 'gbif'
comb_2020 <- comb_2020 %>%
  group_by(powoID) %>%
  filter(n() == 2 &
           all(c("curated", "gbif") %in% workflow)) %>%
  ungroup()

# Reshape to wide format for paired test
wide_df <- comb_2020 %>%
  pivot_wider(names_from = workflow, values_from = cell_count)

# data ready, now for the test

# Paired t-test on log-transformed data
t.test(log(wide_df$gbif), log(wide_df$curated), paired = TRUE)

# Calculate effect size (log ratio)
mean(log(wide_df$gbif) - log(wide_df$curated))

# Or back-transform to get the ratio
exp(mean(log(wide_df$gbif) - log(wide_df$curated)))

# density
p_dens <- ggplot(comb_2020, aes(x = cell_count, fill = workflow, colour = workflow)) +
  geom_density(alpha = 0.5) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Range Size (cells)", y = "Density") +
  scale_fill_manual(values = c("curated" = "#CC79A7", "gbif" = "#509E2F")) +
  scale_color_manual(values = c("curated" = "#CC79A7", "gbif" = "#509E2F")) +
  theme_base(base_size = 14) +
  theme(
    legend.position = "right",
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5)
  )

p_dens

# hex
p_hex <- ggplot(wide_df, aes(x = curated, y = gbif)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "red", linewidth = 0.8) +
  geom_hex(bins = 30, alpha = 0.7) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  scale_fill_viridis_c(name = "Count") +
  labs(x = "Range Size - Curated (no. cells)",
       y = "Range Size - GBIF (no. cells)") +
  theme_base(base_size = 14) +
  theme(
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  )
coord_equal()
getwd()

(p_hex / p_dens)

combo <- cowplot::plot_grid(p_hex, p_dens, labels = c('A', 'B'), nrow = 2)

#ggsave(paste0(getwd(),"/02_figures/Figure_S3_GBIF_vs_Curated.png"), width = 6, height = 5)
ggsave(paste0(getwd(),"/02_figures/Figure_S3_GBIF_vs_Curated.png"), width = 6, height = 8)

