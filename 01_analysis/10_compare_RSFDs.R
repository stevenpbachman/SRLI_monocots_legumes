
# Compare range size frequency distributions

# get the RSFD for all monocots and legumes (WCVP)
# get the RSFD for sampled monocots and legumes
# get the RSFD for Red List Monocots and Legumes

# test a - compare the assessed sample with NE sample
# test b - compare full sample with the Red List
# test c - compare WCVP (truth) with full sample

library(dplyr)
library(readr)

# organise data
# start from scratch so that this is standalone
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# main files
tdwg_level3 <- readRDS(paste0(work_dir, "00_raw_data/tdwg_level3.rds"))
SRLI_dists <- read_csv(paste0(work_dir,"00_raw_data/SRLI_dists.csv"))
distributions <- rWCVPdata::wcvp_distributions
tax_map <- rWCVP::taxonomic_mapping

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

# join names and distributions and filter to the accepted, native etc.
dists <- rWCVPdata::wcvp_names %>% dplyr::select(plant_name_id,powo_id, taxon_name, taxon_status, taxon_rank) %>%
  dplyr::left_join(distributions, by = "plant_name_id") %>%
  dplyr::filter(complete.cases(plant_locality_id),
                taxon_status == "Accepted",
                taxon_rank == "Species",
                location_doubtful == "0",
                extinct == "0",
                introduced == "0")

all_wcvp <- dists %>%
group_by(powo_id) %>%
  count() %>%
  mutate(dataset = "WCVP all")


# Test a.assessed sample vs full sample ----
# Monocots
# 1a. all monocots in the sample
mons_SRLI_assessed <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
   group_by(WCVP_current_id) %>%
   count() %>%
   mutate(dataset = "Assessed SRLI sample")%>%
  mutate(group = "Monocots")

# 1b sampled monocots NE only
mons_SRLI_NE <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 == "NE") %>%
  group_by(WCVP_current_id) %>%
  count() %>%
  mutate(dataset = "NE SRLI sample")%>%
  mutate(group = "Monocots")

# 2. Compare assessed sample with full sample
# Legumes
legs_SRLI_assessed <- SRLI_dists %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
   group_by(WCVP_current_id) %>%
   count() %>%
   mutate(dataset = "Assessed SRLI sample") %>%
  mutate(group = "Legumes")

# 1b sampled monocots NE only
legs_SRLI_NE <- SRLI_dists %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include",
         RL_2020 == "NE") %>%
  group_by(WCVP_current_id) %>%
  count() %>%
  mutate(dataset = "NE SRLI sample") %>%
  mutate(group = "Legumes")

test_a_comb <- bind_rows(mons_SRLI_assessed, mons_SRLI_NE, legs_SRLI_assessed, legs_SRLI_NE)

# plot
# first plot frequency
a <- ggplot(test_a_comb, aes(x = n, fill = interaction(dataset, group), colour = group, linetype = dataset)) +
  geom_density(alpha = 0.3, adjust = 1.5, linewidth = 0.8) +
  scale_x_log10() +
  facet_wrap(~ group) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  scale_fill_manual(values = c(
    "Assessed SRLI sample.Monocots" = "#0072B2",
    "NE SRLI sample.Monocots" = "#0072B2",
    "Assessed SRLI sample.Legumes"  = "#E69F00",
    "NE SRLI sample.Legumes"  = "#E69F00"
  )) +
  scale_linetype_manual(values = c("Assessed SRLI sample" = "solid", "NE SRLI sample" = "dashed")) +
  labs(x = "Number of WGSRPD regions", y = "Density", fill = NULL, colour = NULL, linetype = NULL) +
  guides(
    colour   = "none",
    fill     = "none",
    linetype = guide_legend()
  ) +
  theme_base(base_size = 14) +
  theme(legend.position="bottom",
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ))

ggsave("02_figures/Figure_S7_A.png", width = 6, height = 3)

sum_stat_a <- test_a_comb %>%
  group_by(group,dataset) %>%
  summarise(
    n_species = n(),
    mean_range = mean(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )

# Test b. full sample vs red list ----
# 1. all monocots on Red List
mons_red_list <- read_csv(paste0(work_dir,"00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(className == "LILIOPSIDA") %>%
  left_join(dists, by = c("ipni_id" = "powo_id"), relationship = "many-to-many")

mons_red_list <- mons_red_list %>%
  # Remove Data Deficient
  filter(redlistCategory != "Data Deficient") %>%
  # Create threat_status column
  mutate(threat_status = case_when(
    redlistCategory %in% c("Critically Endangered", "Endangered", "Vulnerable", "Extinct", "Extinct in the Wild") ~ "Threatened",
    TRUE ~ "Not Threatened"
  ))

mons_red_list <- mons_red_list %>%
     group_by(ipni_id, threat_status) %>%
     count() %>%
     mutate(dataset = "Red List") %>%
   filter(complete.cases(ipni_id)) %>%
  mutate(group = "Monocots")


# 2. all legumes on Red List
legs_red_list <- read_csv(paste0(work_dir,"00_raw_data/WCVPv12-RL2023-1.csv")) %>%
  dplyr::filter(familyName == "FABACEAE") %>%
  left_join(dists, by = c("ipni_id" = "powo_id"), relationship = "many-to-many")

legs_red_list <- legs_red_list %>%
  # Remove Data Deficient
  filter(redlistCategory != "Data Deficient") %>%
  # Create threat_status column
  mutate(threat_status = case_when(
    redlistCategory %in% c("Critically Endangered", "Endangered", "Vulnerable", "Extinct", "Extinct in the Wild") ~ "Threatened",
    TRUE ~ "Not Threatened"
  ))

legs_red_list <- legs_red_list %>%
     group_by(ipni_id, threat_status) %>%
     count() %>%
     mutate(dataset = "Red List") %>%
   filter(complete.cases(ipni_id)) %>%
  mutate(group = "Legumes")


combinedRL <- bind_rows(legs_red_list, mons_red_list)
combinedRL <- combinedRL %>%
  ungroup() %>%
  select(c(-ipni_id, -threat_status))

# 3. full sample
mons_SRLI <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include") %>%
   group_by(WCVP_current_id) %>%
   count() %>%
   mutate(dataset = "Full SRLI sample")%>%
  mutate(group = "Monocots")

legs_SRLI <- SRLI_dists %>%
  filter(group == "Legumes",
         SRLI_2024_analysis == "include") %>%
   group_by(WCVP_current_id) %>%
   count() %>%
   mutate(dataset = "Full SRLI sample") %>%
  mutate(group = "Legumes")

full_sample <- bind_rows(mons_SRLI,legs_SRLI)
full_sample <- full_sample %>%
  ungroup() %>%
  select(-WCVP_current_id)

test_b_comb <- bind_rows(full_sample, combinedRL)

# plot
ggplot(test_b_comb, aes(x = n, fill = interaction(dataset, group), colour = group, linetype = dataset)) +
  geom_density(alpha = 0.3, adjust = 1.5, linewidth = 0.8) +
  scale_x_log10() +
  facet_wrap(~ group) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  scale_fill_manual(values = c(
    "Red List.Monocots" = "#0072B2",
    "Full SRLI sample.Monocots" = "#0072B2",
    "Red List.Legumes"  = "#E69F00",
    "Full SRLI sample.Legumes"  = "#E69F00"
  )) +
  scale_linetype_manual(values = c("Red List" = "solid", "Full SRLI sample" = "dashed")) +
  labs(x = "Number of WGSRPD regions", y = "Density", fill = NULL, colour = NULL, linetype = NULL) +
  guides(
    colour   = "none",
    fill     = "none",
    linetype = guide_legend()
  ) +
  theme_base(base_size = 14) +
  theme(legend.position="bottom",
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ))

#ggsave("02_figures/Figure_S7_B.png", width = 6, height = 3)

sum_stat_b <- test_b_comb %>%
  group_by(group,dataset) %>%
  summarise(
    n_species = n(),
    #median_range = median(n),
    mean_range = mean(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )


# Test c. all WCVP vs full sample and red list ----
# MONOCOTS -
# 1. all monocots from WCVP -
# get list of monocot families
APGIV_families_orders_IUCN_crosswalk <- read.csv(paste0(work_dir,"00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))
apg_mon_fams <- APGIV_families_orders_IUCN_crosswalk %>%
  filter(HigherGroups == "Monocots")

wcvp_mons <- wcvp_acc %>%
  left_join(apg_mon_fams, by = c("family" = "APG.IV.Family")) %>%
  filter(HigherGroups == "Monocots")

# all WCVP monocots geography - RSFD
mons_geog_wcvp <- wcvp_mons %>%
  left_join(dists, by = "powo_id") %>%
  group_by(powo_id) %>%
  count()

mons_WCVP <- mons_geog_wcvp %>%
  mutate(dataset = "WCVP") %>%
  mutate(group = "Monocots")

#  Legumes -
# 1. all legumes from WCVP -
wcvp_legs <- wcvp_acc %>%
  filter(family == "Fabaceae")

# all WCVP legumes geography - RSFD
legs_geog_wcvp <- wcvp_legs %>%
  left_join(dists, by = "powo_id") %>%
  group_by(powo_id) %>%
  count()

legs_WCVP <- legs_geog_wcvp %>%
  mutate(dataset = "WCVP") %>%
  mutate(group = "Legumes")

wcvp <- bind_rows(legs_WCVP, mons_WCVP)

test_c_comb <- bind_rows(wcvp, combinedRL, full_sample)

# combined plot
# plot
# plot
c <- ggplot(test_c_comb, aes(x = n, fill = group, colour = group)) +
  geom_density(alpha = 0.3, adjust = 1.5, linewidth = 0.8) +
  scale_x_log10() +
  facet_grid(group ~ dataset) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  scale_fill_manual(values  = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  labs(x = "Number of WGSRPD regions", y = "Density", fill = NULL, colour = NULL) +
  guides(colour = "none", fill = "none") +
  theme_base(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5)
  )

ggsave("02_figures/Figure_S7_c.png", width = 6, height = 4)

sum_stat_c <- test_c_comb %>%
  group_by(group,dataset) %>%
  summarise(
    n_species = n(),
    #median_range = median(n),
    mean_range = mean(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )

# combined plot Figure_S7_rsfds.png
combo <- cowplot::plot_grid(a, c, labels = c('A', 'B'), nrow = 2)
ggsave(("02_figures/Figure_S7_rsfds.png"),width = 21, height = 20,
       units = "cm")

# combined NE sample vs assessed sample
# Red list vs full sample vs WCVP
# Table S10 ----
# save Table_S10_rsfds.csv




# year of description ----
# get WCVP list of accepted species and distributions ----

# filter to accepted species
df <- wcvp_names %>%
  filter(taxon_status=="Accepted",
         taxon_rank == "Species",
         is.na(genus_hybrid),
         is.na(species_hybrid),
         is.na(infraspecies)) #should be redundant but better to be safe

df_nobasio <- df %>% select(plant_name_id, taxon_name, family, genus,
                            lifeform_description, first_published,
                            basionym_plant_name_id) %>%
  filter(is.na(basionym_plant_name_id))


df_basio <- df %>% select(plant_name_id, taxon_name, family, genus,
                          lifeform_description, first_published,
                          basionym_plant_name_id) %>%
  filter(!is.na(basionym_plant_name_id)) %>%
  left_join(rWCVPdata::wcvp_names %>%
              select(plant_name_id, basio_year=first_published),
            by=c("basionym_plant_name_id"="plant_name_id"))

df_basio$year <- df_basio$basio_year
df_nobasio$year <- df_nobasio$first_published

df <- bind_rows(df_nobasio, df_basio)

df$year <- gsub("\\(", "",df$year )
df$year <- gsub("\\)", "",df$year )
years <- str_extract_all(df$year, "\\d\\d\\d\\d")
df$year <- unlist(lapply(years, function(x) x[1]))
df$year <- as.numeric(df$year)

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
         is.na(species_hybrid)) #|>
  #select(powo_id, taxon_name, family, genus, lifeform_description, climate_description)
#nrow(wcvp_acc)

# join names and distributions and filter to the accepted, native etc.
dists <- rWCVPdata::wcvp_names %>% dplyr::select(plant_name_id,powo_id, taxon_name, taxon_status, taxon_rank) %>%
  dplyr::left_join(distributions, by = "plant_name_id") %>%
  dplyr::filter(complete.cases(plant_locality_id),
                taxon_status == "Accepted",
                taxon_rank == "Species",
                location_doubtful == "0",
                extinct == "0",
                introduced == "0")

all_wcvp <- dists %>%
group_by(powo_id, plant_name_id) %>%
#  group_by(powo_id) %>%
  count() %>%
  mutate(dataset = "WCVP all")

# join ids to distributions
#all_wcvp to df
wcvp_years <- all_wcvp %>% left_join(df, by = "plant_name_id")

# filter on monocots and legumes
# 1. all monocots from WCVP -
# get list of monocot families
APGIV_families_orders_IUCN_crosswalk <- read.csv(paste0(work_dir,"00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))
apg_mon_fams <- APGIV_families_orders_IUCN_crosswalk %>%
  filter(HigherGroups == "Monocots")

wcvp_years_mons <- wcvp_years %>%
  left_join(apg_mon_fams, by = c("family" = "APG.IV.Family")) %>%
  filter(HigherGroups == "Monocots")%>%
  mutate(group = "Monocots")

# legumes
wcvp_years_legs <- wcvp_years %>%
  filter(family == "Fabaceae") %>%
  mutate(group = "Legumes")

comb_years <- bind_rows(wcvp_years_mons, wcvp_years_legs) %>%
  mutate(year = if_else(year < 2006, "pre-2006", "post-2006")) %>%
  dplyr::filter(complete.cases(year))

#tes <- comb_years %>% group_by(year) %>% count()

#plot
c <- ggplot(comb_years, aes(x = n, fill = interaction(year, group), colour = group, linetype = year)) +
  geom_density(alpha = 0.3, adjust = 1.5, linewidth = 0.8) +
  scale_x_log10() +
  facet_wrap(year ~ group) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  scale_fill_manual(values = c(
    "pre-2006.Monocots" = "#0072B2",
    "post-2006.Monocots" = "#0072B2",
    "pre-2006.Legumes"  = "#E69F00",
    "post-2006.Legumes"  = "#E69F00"
  )) +
  scale_linetype_manual(values = c("Red List" = "solid", "Full SRLI sample" = "dashed")) +
  labs(x = "Number of WGSRPD regions", y = "Density", fill = NULL, colour = NULL, linetype = NULL) +
  guides(
    colour   = "none",
    fill     = "none",
    linetype = guide_legend()
  ) +
  theme_base(base_size = 14) +
  theme(legend.position="bottom",
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ))
c

sum_stat_yod <- comb_years %>%
  group_by(group,year) %>%
  summarise(
    n_species = n(),
    #median_range = median(n),
    mean_range = mean(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )
sum_stat_yod$dataset <- sum_stat_yod$year
sum_stat_yod <- subset(sum_stat_yod, select = -year)

# combined
Table_S10_rsfds <- bind_rows(sum_stat_a, sum_stat_c, sum_stat_yod)
write.csv(Table_S10_rsfds,
          "02_figures/Table_S10_rsfds.csv",
          row.names = FALSE)
