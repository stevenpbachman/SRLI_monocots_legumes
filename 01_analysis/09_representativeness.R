
# Representativeness ----

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
library(ggpubr)
library(hexbin)
library(ggthemes)
library(stringr)
library(purrr)
library(rCAT)

# Figure S5 - representativeness of sample ----
# add validation here? representativeness of GBIF validation sample

# organise data
# start from scratch so that this is standalone
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir, "00_raw_data/SRLI_2026.csv")) %>%
  dplyr::filter(SRLI_2024_analysis == "include")

# 1 - taxonomic ----
# 1a - taxonomic - full original sample
# Load in the wcvp
wcvp_names <- read_delim(paste0(work_dir,"00_raw_data/wcvp_names.csv"),
                         delim = "|",
                         escape_double = FALSE,
                         trim_ws = TRUE)

# get the number of accepted species
wcvp_acc <-
  wcvp_names |>
  filter(taxon_rank == "Species",
         taxon_status == "Accepted",
         !is.na(taxon_name),
         is.na(genus_hybrid),
         is.na(species_hybrid)) |>
  select(powo_id, taxon_name, family, genus, lifeform_description, climate_description)
nrow(wcvp_acc)

###  load SRLI names list
#SRLI_2024 <- read_csv("00_raw_data/SRLI_2024.csv")

comb_val <- read_csv("02_figures/Table_S7_validation_tests.csv")

# validation data - at 1km
validated_sp <- comb_val %>%
  dplyr::filter(method == "1km", valid_test == "TRUE") %>%
  group_by(sp_id) %>%
  count() %>%
  left_join(wcvp_acc, by = c("sp_id" = "powo_id")) %>%
  ungroup() %>%
  select(-sp_id)

##### WCVP legumes
# extract accepted fabaceae only
Fabaceae <- wcvp_acc %>%
  dplyr::filter(family == "Fabaceae")

# get the species count per genera in Fabaceae
wcvp_fab_genera = Fabaceae %>%
  count(genus) %>%
  arrange(desc(n))

##### SRLI legumes
srli_fabs_included <- SRLI_2026 %>%
  dplyr::filter(family == "Fabaceae",
                SRLI_2024_analysis == "include") %>%
  count(genus) %>%
  arrange(desc(n)) %>%
  mutate(source = "included_sample")

srli_fabs_assessed <- SRLI_2026 %>%
  filter(family == "Fabaceae",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  count(genus) %>%
  arrange(desc(n)) %>%
  mutate(source = "assessed_sample")

red_list_fabs <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(familyName == "FABACEAE") %>%
  count(genusName) %>%
  arrange(desc(n)) %>%
  mutate(source = "red_list") %>%
  rename("genus" = "genusName" )

validated_fabs <- validated_sp %>%
  dplyr::filter(family == "Fabaceae") %>%
  count(genus) %>%
  arrange(desc(n)) %>%
  mutate(source = "validated")

# merge them with WCVP
fabs_merged <- bind_rows(srli_fabs_included, srli_fabs_assessed, red_list_fabs, validated_fabs)

# merge the genera tables
fab_gen_merge <- inner_join(fabs_merged, wcvp_fab_genera, by = "genus")

# factors to set order
fab_gen_merge <- fab_gen_merge %>%
  mutate(source = factor(source, levels = c("included_sample","assessed_sample", "validated","red_list")))

# plot
legs <- ggplot(fab_gen_merge, aes(x = n.x, y = n.y, group = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), # Map shape to source
             color = "grey22", stroke = 0.2) +
  geom_smooth(aes(color = source), method = "lm", se = TRUE) + # Ensure smooth lines match colors
  scale_fill_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                               "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                                "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21,  # Circle with fill
                                "assessed_sample" = 22,  # Square with fill
                                "included_sample" = 24)) + # Triangle with fill
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 2800, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Species per Legume genera in samples",
       y = "Species per Legume genera in WCVP")
legs

#ggsave(("02_figures/Figure_S1_legumes.png"),width = 20, height = 9, units = "cm")

##### WCVP monocots
# get list of monocot families
APGIV_families_orders_IUCN_crosswalk <- read.csv(paste0(work_dir, "00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))
apg_mon_fams <- APGIV_families_orders_IUCN_crosswalk %>%
  filter(HigherGroups == "Monocots")

wcvp_mons <- wcvp_acc %>%
  left_join(apg_mon_fams, by = c("family" = "APG.IV.Family")) %>%
  filter(HigherGroups == "Monocots")

wcvp_mons_fams <- wcvp_mons |>
  group_by(family) %>%
  count()

# get the count of monocot families in the sample
srli_mons_included = SRLI_2026 %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include") %>%
  count(family) %>%
  arrange(desc(n)) %>%
  mutate(source = "included_sample")

srli_mons_assessed <- SRLI_2026 %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  count(family) %>%
  arrange(desc(n)) %>%
  mutate(source = "assessed_sample")

red_list_mons <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(className == "LILIOPSIDA") %>%
  count(familyName) %>%
  arrange(desc(n)) %>%
  mutate(source = "red_list") %>%
  rename("family" = "familyName" )

red_list_mons$family <- str_to_lower(red_list_mons$family)
red_list_mons$family <- str_to_title(red_list_mons$family)

validated_mons <- validated_sp %>%
  dplyr::filter(family != "Fabaceae") %>%
  group_by(family) %>%
  count() %>%
  mutate(source = "validated")

# merge them with WCVP
mons_merged <- bind_rows(srli_mons_included, srli_mons_assessed, validated_mons, red_list_mons)

# merge the genera tables
mons_fam_merge <- inner_join(mons_merged, wcvp_mons_fams, by = "family")

# factors to set order
mons_fam_merge <- mons_fam_merge %>%
  mutate(source = factor(source, levels = c("included_sample","assessed_sample", "validated", "red_list")))

# plot
mons <- ggplot(mons_fam_merge, aes(x = n.x, y = n.y, group = source, color = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), color = "grey22", stroke = 0.2) +
  geom_smooth(aes(color = source), method = "lm", se = TRUE) + # Ensure smooth lines match colors
  scale_fill_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                               "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                                "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21, "assessed_sample" = 22,  "included_sample" = 24)) + # Triangle with fill
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 28000, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Species per Monocot families in samples",
       y = "Species per Monocot families in WCVP")
mons

#ggsave(("02_figures/Figure_S1_mons.png"),width = 20, height = 9, units = "cm")

combo <- cowplot::plot_grid(legs, mons, labels = c('A', 'B'), nrow = 2)

ggsave(("02_figures/Figure_S5_AB_tax.png"),width = 21, height = 20, units = "cm")


# 2 - lifeform ----
# pull data together to make long
# life forms per sample: i) original sample, ii) assessed sample iii) red list sample
# monocots

lifeform_mapping_2023_matilda_update <- read.csv(paste0(work_dir, "00_raw_data/lifeform_mapping_2023_matilda_update.csv"))
lf <- lifeform_mapping_2023_matilda_update

APGIV_families_orders_IUCN_crosswalk <- read.csv(paste0(work_dir, "00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))
apg_mon_fams <- APGIV_families_orders_IUCN_crosswalk %>%
  filter(HigherGroups == "Monocots")

wcvp_acc <-
  wcvp_names |>
  filter(taxon_rank == "Species",
         taxon_status == "Accepted",
         !is.na(taxon_name),
         is.na(genus_hybrid),
         is.na(species_hybrid)) |>
  select(powo_id, taxon_name, family, genus, lifeform_description, climate_description)

wcvp_mons <- wcvp_acc %>%
  left_join(apg_mon_fams, by = c("family" = "APG.IV.Family")) %>%
  filter(HigherGroups == "Monocots")

lifeform <- left_join(wcvp_mons, lf, by="lifeform_description", relationship = "many-to-many") |>
  select(powo_id, latest_lifeform) |>
  filter(complete.cases(latest_lifeform))

# all life forms in the world checklist -
wcvp_mons_lf <- lifeform |>
  group_by(latest_lifeform) %>%
  count()

validated_sp <- comb_val %>%
  dplyr::filter(method == "1km", valid_test == "TRUE") %>%
  group_by(sp_id) %>%
  count() %>%
  left_join(wcvp_acc, by = c("sp_id" = "powo_id"))

# now get the samples
# life form from original sample, assessed sample and red list
SRLI_2026_lf_included <- SRLI_2026 %>%
  left_join(lifeform, by = c("WCVP_current_id" = "powo_id"), relationship = "many-to-many") %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "included_sample")

SRLI_2026_lf_assessed <- SRLI_2026 %>%
  left_join(lifeform, by = c("WCVP_current_id" = "powo_id"), relationship = "many-to-many") %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "assessed_sample")

SRLI_2026_lf_red_list <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(className == "LILIOPSIDA") %>%
  left_join(lifeform, by = c("ipni_id" = "powo_id"), relationship = "many-to-many") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "red_list")

SRLI_2026_lf_validated <- validated_sp %>%
  left_join(lifeform, by = c("sp_id" = "powo_id"), relationship = "many-to-many") %>%
  dplyr::filter(family != "Fabaceae") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "validated")

# merge them with WCVP
mons_lf_merged <- bind_rows(SRLI_2026_lf_included, SRLI_2026_lf_assessed, SRLI_2026_lf_validated, SRLI_2026_lf_red_list)

# merge the genera tables
mons_lf_merge <- inner_join(mons_lf_merged, wcvp_mons_lf, by = "latest_lifeform")

# factors to set order
mons_lf_merge <- mons_lf_merge %>%
  mutate(source = factor(source, levels = c("included_sample","assessed_sample", "validated", "red_list")))

# plot
mons_lf <- ggplot(mons_lf_merge, aes(x = n.x, y = n.y, group = source, color = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), # Map shape to source
             color = "grey22", stroke = 0.2) +
  scale_fill_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                               "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                                "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21, "assessed_sample" = 22,  "included_sample" = 24)) + # Triangle with fill
  geom_smooth(method = "lm", se=TRUE) +
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 55000, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Monocot plant life form in samples",
       y = "Monocot plant life form in WCVP")
mons_lf

#ggsave(("02_figures/Extended_Data_Figure_S1_monocots_lifeform.png"),width = 20, height = 15, units = "cm")

# legumes
lifeform_fab <- left_join(Fabaceae, lf, by="lifeform_description", relationship = "many-to-many") |>
  select(powo_id, latest_lifeform) |>
  filter(complete.cases(latest_lifeform))

# all life forms in the world checklist -
wcvp_fabs_lf <- lifeform_fab |>
  group_by(latest_lifeform) %>%
  count()

# now get the samples
# life form from original sample, assessed sample and red list
SRLI_2026_lf_included <- SRLI_2026 %>%
  left_join(lifeform_fab, by = c("WCVP_current_id" = "powo_id"), relationship = "many-to-many") %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "included_sample")

SRLI_2026_lf_assessed <- SRLI_2026 %>%
  left_join(lifeform_fab, by = c("WCVP_current_id" = "powo_id"), relationship = "many-to-many") %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "assessed_sample")

SRLI_2026_lf_red_list <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(familyName == "FABACEAE") %>%
  left_join(lifeform_fab, by = c("ipni_id" = "powo_id"), relationship = "many-to-many") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "red_list")

SRLI_2026_lf_validated <- validated_sp %>%
  left_join(lifeform_fab, by = c("sp_id" = "powo_id"), relationship = "many-to-many") %>%
  dplyr::filter(family == "Fabaceae") %>%
  group_by(latest_lifeform) %>%
  count() %>%
  mutate(source = "validated")

# merge them with WCVP
fabs_lf_merged <- bind_rows(SRLI_2026_lf_included, SRLI_2026_lf_assessed, SRLI_2026_lf_validated, SRLI_2026_lf_red_list)

# merge the genera tables
fabs_lf_merge <- inner_join(fabs_lf_merged, wcvp_fabs_lf, by = "latest_lifeform")

# factors to set order
fabs_lf_merge <- fabs_lf_merge %>%
  mutate(source = factor(source, levels = c("included_sample","assessed_sample", "validated", "red_list")))

# plot
fabs_lf <- ggplot(fabs_lf_merge, aes(x = n.x, y = n.y, group = source, color = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), # Map shape to source
             color = "grey22", stroke = 0.2) +
  geom_smooth(aes(color = source), method = "lm", se = TRUE) + # Ensure smooth lines match colors
  scale_fill_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                               "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                                "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21,  # Circle with fill
                                "assessed_sample" = 22,  # Square with fill
                                "included_sample" = 24)) + # Triangle with fill
  geom_smooth(method = "lm", se=TRUE) +
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 50000, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Legumes plant life form in samples", #http://127.0.0.1:47533/graphics/plot_zoom_png?width=617&height=189
       y = "Legumes plant life form in WCVP")
fabs_lf

combo_lf <- cowplot::plot_grid(fabs_lf, mons_lf, labels = c('C', 'D'), nrow = 2)
ggsave(("02_figures/Figure_S5_CD_lifeform.png"),width = 20, height = 20, units = "cm")


# 3 - geographic ----
tdwg_level3 <- readRDS(paste0(work_dir, "00_raw_data/tdwg_level3.rds"))
SRLI_dists <- read_csv(paste0(work_dir, "00_raw_data/SRLI_dists.csv"))
distributions <- rWCVPdata::wcvp_distributions

# join names and distributions and filter to the accepted, native etc.
dists <- wcvp_names %>% dplyr::select(plant_name_id,powo_id, taxon_name, taxon_status, taxon_rank) %>%
  dplyr::left_join(distributions, by = "plant_name_id") %>%
  dplyr::filter(complete.cases(plant_locality_id),
                taxon_status == "Accepted",
                taxon_rank == "Species",
                location_doubtful == "0",
                extinct == "0",
                introduced == "0")

# all monocots included - sample
mons_geog_included <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "included_sample")

# all monocots assessed - sample
mons_geog_assesed <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "assessed_sample")

# monocot RL geog
mons_geog_wcvp <- wcvp_mons %>%
  left_join(dists, by = "powo_id")

# all monocots red list
SRLI_2024_geog_red_list <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  left_join(mons_geog_wcvp, by = c("ipni_id" = "powo_id"), relationship = "many-to-many") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "red_list")

# monocot validated geog
validated_geog <- validated_sp %>%
  left_join(dists, by = c("sp_id" = "powo_id")) %>%
  dplyr::filter(family != "Fabaceae") %>%
  select(area_code_l3, sp_id) %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "validated")

# all WCVP monocots
mons_geog_wcvp <- wcvp_mons %>%
  left_join(dists, by = "powo_id") %>%
  group_by(area_code_l3) %>%
  count()

# bind assessed and Red list
mons_geog_merged <- bind_rows(mons_geog_included, mons_geog_assesed, validated_geog, SRLI_2024_geog_red_list)

# join with WCVP geog
mons_geog_merged <- inner_join(mons_geog_merged, mons_geog_wcvp, by = "area_code_l3")

# factors to set order
mons_geog_merged <- mons_geog_merged %>%
  mutate(source = factor(source, levels = c("included_sample", "assessed_sample", "validated", "red_list"))) %>%
  filter(complete.cases(area_code_l3))

# plot
mons_geog <- ggplot(mons_geog_merged, aes(x = n.x, y = n.y, group = source, color = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), color = "grey22", stroke = 0.2) +
  geom_smooth(aes(color = source), method = "lm", se = TRUE) + # Ensure smooth lines match colors
  scale_fill_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                               "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#87BFDA", "red_list" = "#E6F3FA",
                                "assessed_sample" = "#4297C3","included_sample" = "#0072B2")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21, "assessed_sample" = 22,  "included_sample" = 24)) + # Triangle with fill
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 7500, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Monocot species per L3 region in samples",
       y = "Monocot species per L3 region in WCVP")
mons_geog

#ggsave(("02_figures/Extended_Data_Figure_S1_monocots_geography.png"),width = 20, height = 15, units = "cm")



# legumes - geog

# all legumes included - sample
legs_geog_included <- SRLI_dists %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "included_sample")

# all legumes assessed - sample
legs_geog_assesed <- SRLI_dists %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "assessed_sample")

# legumes RL geog
legs_geog_wcvp <- Fabaceae %>%
  left_join(dists, by = "powo_id")

SRLI_2024_geog_red_list <- read_csv(paste0(work_dir, "00_raw_data/WCVPv12-RL2023-1.csv"))%>%
  left_join(legs_geog_wcvp, by = c("ipni_id" = "powo_id"), relationship = "many-to-many") %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "red_list")

# legumes validated geog
validated_geog <- validated_sp %>%
  left_join(dists, by = c("sp_id" = "powo_id")) %>%
  dplyr::filter(family == "Fabaceae") %>%
  select(area_code_l3, sp_id) %>%
  group_by(area_code_l3) %>%
  count() %>%
  mutate(source = "validated")

# all WCVP legumes
legs_geog_wcvp <- Fabaceae %>%
  left_join(dists, by = "powo_id") %>%
  group_by(area_code_l3) %>%
  count()

# merge them with WCVP
legs_geog_merged <- bind_rows(legs_geog_included, legs_geog_assesed, validated_geog, SRLI_2024_geog_red_list)

# merge the tables
legs_geog_merged <- inner_join(legs_geog_merged, legs_geog_wcvp, by = "area_code_l3")

# factors to set order
legs_geog_merged <- legs_geog_merged %>%
  mutate(source = factor(source, levels = c("included_sample", "assessed_sample", "validated", "red_list"))) %>%
  filter(complete.cases(area_code_l3))

# plot
legs_geog <- ggplot(legs_geog_merged, aes(x = n.x, y = n.y, group = source, color = source)) +
  geom_point(aes(size = 3, fill = source, shape = source), # Map shape to source
             color = "grey22", stroke = 0.2) +
  geom_smooth(aes(color = source), method = "lm", se = TRUE) + # Ensure smooth lines match colors
  scale_fill_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                               "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Fill scale for points
  scale_color_manual(values = c("validated" = "#F8C98E", "red_list" = "#FFE5CC",
                                "assessed_sample" = "#F0B347","included_sample" = "#E69F00")) + # Color scale for smooth lines
  scale_shape_manual(values = c("validated" = 23, "red_list" = 21,  # Circle with fill
                                "assessed_sample" = 22,  # Square with fill
                                "included_sample" = 24)) + # Triangle with fill
  theme_base(base_size = 12) +
  facet_wrap(vars(source), scales = "free_x", nrow = 1) +
  theme(legend.position="none") +
  stat_regline_equation(label.y = 1500, aes(label = ..rr.label..), size = 5, label.x.npc = "left", color = "black") +
  labs(x = "Legume species per L3 region in samples",
       y = "Legume species per L3 region in WCVP")
legs_geog

combo <- cowplot::plot_grid(legs_geog, mons_geog, labels = c('E', 'F'), nrow = 2)

ggsave(("02_figures/Figure_S5_EF_geog.png"),width = 20, height = 15, units = "cm")


# new test Chi ----
# For one source subset at a time
df_sub <- fab_gen_merge |> filter(source == "included_sample")

total_sample <- sum(df_sub$n.x)
total_full   <- sum(df_sub$n.y)

# Expected counts in the sample if proportions matched full set
expected <- total_sample * (df_sub$n.y / total_full)

chisq.test(x = df_sub$n.x, p = df_sub$n.y / total_full)


test <- chisq.test(x = df_sub$n.x, p = df_sub$n.y / sum(df_sub$n.y))

# Standardised residuals > |2| are the problem genera
residuals <- data.frame(
  genus    = df_sub$genus,
  observed = df_sub$n.x,
  expected = test$expected,
  residual = test$stdres
) |>
  arrange(desc(abs(residual)))

head(residuals, 20)
