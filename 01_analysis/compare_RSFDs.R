
# Compare range size frequency distributions

# get the RSFD for all monocots and legumes (WCVP)
# get the RSFD for sampled monocots and legumes
# get the RSFD for Red List Monocots and Legumes

# compare sampled and red list with all (WCVP)

# main files
tdwg_level3 <- readRDS("00_raw_data/tdwg_level3.rds")
SRLI_dists <- read_csv("00_raw_data/SRLI_dists.csv")
distributions <- rWCVPdata::wcvp_distributions
tax_map <- rWCVP::taxonomic_mapping

# Load in the wcvp
wcvp_names <- read_delim("00_raw_data/wcvp_names.csv",
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
  mutate(dataset = "WCVP_all")

# MONOCOTS ----
# 1. all monocots from WCVP ----
# get list of monocot families
APGIV_families_orders_IUCN_crosswalk <- read.csv("00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv")
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
  mutate(dataset = "WCVP_mons")

# 2. all monocots on Red List ----
mons_red_list <- read_csv("00_raw_data/WCVPv12-RL2023-1.csv") %>%
  dplyr::filter(className == "LILIOPSIDA") %>%
  left_join(dists, by = c("ipni_id" = "powo_id"), relationship = "many-to-many")
#%>%

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
     mutate(dataset = "red_list") %>%
   filter(complete.cases(ipni_id))

# 3. all monocots in the sample ----
mons_SRLI_assessed <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 != "NE") %>%
   group_by(WCVP_current_id) %>%
   count() %>%
   mutate(dataset = "assessed_sample")

# 3b sampled monocots NE only
mons_SRLI_NE <- SRLI_dists %>%
  filter(group == "Monocots",
         SRLI_2024_analysis == "include",
         RL_2020 == "NE") %>%
  group_by(WCVP_current_id) %>%
  count() %>%
  mutate(dataset = "NE_sample")

srli_comb <- bind_rows(mons_SRLI_assessed, mons_SRLI_NE)

# LEGUMES ----
# do the same for legumes


#combine and plot ----
combined <- bind_rows(mons_red_list, mons_WCVP, mons_SRLI, all_wcvp)

# first plot frequency
ggplot(srli_comb, aes(x = n, fill = dataset, colour = dataset)) +
  geom_density(alpha = 0.3, adjust = 1.5) +
  scale_x_log10() +
  labs(x = "Number of countries", y = "Density", fill = NULL, colour = NULL) +
  theme_minimal()

ggplot(srli_comb, aes(x = n, colour = dataset)) +
  stat_ecdf() +
  scale_x_log10() +
  labs(x = "Number of countries", y = "Cumulative proportion") +
  theme_minimal()

srli_comb %>%
  group_by(dataset) %>%
  summarise(
    n_species = n(),
    median_range = median(n),
    mean_range = mean(n),
    #skewness = moments::skewness(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )

ks.test(mons_red_list$n, mons_WCVP$n)
ks.test(mons_SRLI$n, mons_WCVP$n)
