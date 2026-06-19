
# change in new species

# influence of new species of monocots and legumes on the sample

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

#all_wcvp$plant_name_id <- as.character(all_wcvp$plant_name_id)

# year of description
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

tes <- comb_years %>% group_by(year) %>% count()

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

sum_stat <- comb_years %>%
  group_by(group,year) %>%
  summarise(
    n_species = n(),
    #median_range = median(n),
    mean_range = mean(n),
    pct_singletons = mean(n == 1)  # proportion in just 1 country
  )
sum_stat
