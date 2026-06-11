
# Data coverage summary ----

# summarise data so it can be used for a Sankey?
# first section of results

library(readr)
library(dplyr)
library(tidyverse)
library(networkD3)
library(htmlwidgets)
library(webshot2)

#library(sankeyD3)

work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# Functions ----

# --- Colour helpers ---
node_colour_js <- function(colours) {
  paste0(
    'd3.scaleOrdinal().range(["',
    paste(colours, collapse = '","'),
    '"])'
  )
}

make_node_colours <- function(base_hex, n_stages = 5) {
  # Dark to light tints across stages
  tints <- colorRampPalette(c(base_hex, "#FFFFFF"))(n_stages + 2)[1:n_stages]
  tints
}

flow_to_sankey <- function(flow_data) {

  original     <- sum(flow_data$n)
  accepted     <- sum(flow_data$n[flow_data$accepted != "not_accepted"])
  not_accepted <- sum(flow_data$n[flow_data$accepted == "not_accepted"])
  assessed     <- sum(flow_data$n[flow_data$assessed == "assessed"])
  not_assessed <- sum(flow_data$n[flow_data$assessed == "not_assessed"])
  mapped       <- sum(flow_data$n[flow_data$mapped == "mapped"])
  not_mapped   <- sum(flow_data$n[flow_data$mapped == "not_mapped"])
  curated <- sum(flow_data$n[flow_data$workflow == "curated"], na.rm = TRUE)
  GBIF    <- sum(flow_data$n[flow_data$workflow == "GBIF"], na.rm = TRUE)

  nodes <- data.frame(name = c(
    paste0("Sample (", original, ")"),        # 0
    paste0("Accepted (", accepted, ")"),         # 1
    paste0("Not accepted (", not_accepted, ")"), # 2
    paste0("Assessed (", assessed, ")"),         # 3
    paste0("Not assessed (", not_assessed, ")"), # 4
    paste0("Mapped (", mapped, ")"),             # 5
    paste0("Not mapped (", not_mapped, ")"),      # 6
    paste0("Curated (", curated, ")"),        # 7
    paste0("GBIF (", GBIF, ")")               # 8
    ))

  links <- data.frame(
    source = c(0, 0, 1, 1, 3, 3, 5, 5),
    target = c(2, 1, 4, 3, 6, 5, 7, 8),   # dropped first → peels off top
    value  = c(not_accepted, accepted, not_assessed, assessed, not_mapped, mapped,
               curated, GBIF)
  )

  list(nodes = nodes, links = links)
}

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir,"00_raw_data/SRLI_2026.csv"))

SRLI_2026 %>% group_by(group, SRLI_2024_analysis) %>% count()

alluvial_data <- SRLI_2026 %>% select(group, SRLI_2024_analysis,
                                     SRLI_2020_final_cat, aoh_map,
                                     AOH_map_valid, WCVP_current_id, final_workflow) %>%
  dplyr::filter(SRLI_2024_analysis == "include")

alluvial_data$original <- TRUE

alluvial_data <- alluvial_data %>%
  mutate(accepted = case_when(
    SRLI_2024_analysis == "include" ~ TRUE,
    SRLI_2024_analysis == "exclude" ~ FALSE,
  )) %>%
  mutate(assessed = case_when(
    SRLI_2024_analysis == "include" & SRLI_2020_final_cat != "NE" ~ TRUE,
    SRLI_2024_analysis == "include" & SRLI_2020_final_cat == "NA" ~ FALSE,
    SRLI_2024_analysis == "exclude" ~ NA,
  )) %>%
  mutate(mapped = case_when(
    SRLI_2024_analysis == "include" & aoh_map == "yes" ~ TRUE,
    SRLI_2024_analysis == "include" & aoh_map == "NA" ~ FALSE,
    SRLI_2024_analysis == "exclude" ~ NA,
  )) %>%
  select(original, accepted, assessed, mapped, group, final_workflow)

alluvial_data <- alluvial_data %>%
  mutate(
    workflow = case_when(
      final_workflow == "curated" ~ "curated",
      final_workflow == "GBIF" ~ "GBIF",
      TRUE ~ NA_character_
    )
  )

# split into mono and legs
mono_sank <- alluvial_data %>% filter(group == "Monocots")
legs_sank <- alluvial_data %>% filter(group == "Legumes")

# summary tables for supplementary
comb <- bind_rows(mono_sank,legs_sank)

# data_coverage_summary <-
#   summary_table <- comb %>%
#   group_by(group) %>%
#   summarise(
#     across(
#       where(is.logical),
#       ~ sum(.x, na.rm = TRUE),
#       .names = "{.col}"
#     )
#   )

data_coverage_summary <- comb %>%
  group_by(group) %>%
  summarise(
    across(
      where(is.logical),
      ~ sum(.x, na.rm = TRUE)
    ),
    curated = sum(final_workflow == "curated", na.rm = TRUE),
    GBIF    = sum(final_workflow == "GBIF", na.rm = TRUE)
  )


# write this table to results and merge as a single excel file in sup mat
write.csv(data_coverage_summary,
          "02_figures/Table_S1_data_coverage_summary.csv",
          row.names = FALSE)

# into groups
mon_plot_df <- mono_sank %>%
  mutate(
    original = "Sample",
    accepted = ifelse(accepted, "Accepted", "not_accepted"),
    assessed = case_when(
      accepted == "not_accepted" ~ "not_assessed",
      assessed ~ "assessed",
      TRUE ~ "not_assessed"
    ),
    mapped = case_when(
      assessed == "not_assessed" ~ "not_mapped",
      mapped ~ "mapped",
      TRUE ~ "not_mapped"
    ),
    workflow = case_when(
      mapped == "not_mapped" ~ NA_character_,
      workflow == "curated" ~ "curated",
      workflow == "GBIF" ~ "GBIF",
      TRUE ~ NA_character_
    )
  )

legs_plot_df <- legs_sank %>%
  mutate(
    original = "Sample",
    accepted = ifelse(accepted, "Accepted", "not_accepted"),
    assessed = case_when(
      accepted == "not_accepted" ~ "not_assessed",
      assessed ~ "assessed",
      TRUE ~ "not_assessed"
    ),
    mapped = case_when(
      assessed == "not_assessed" ~ "not_mapped",
      mapped ~ "mapped",
      TRUE ~ "not_mapped"
    ),
    workflow = case_when(
      mapped == "not_mapped" ~ NA_character_,
      workflow == "curated" ~ "curated",
      workflow == "GBIF" ~ "GBIF",
      TRUE ~ NA_character_
    )
  )


# ── MONOCOTS ────────────────────────────────────────────────────────────────
mono_tints <- make_node_colours("#0072B2", 5)

mono_colours <- c(
  mono_tints[1],
  mono_tints[2],
  #dropped_col,
  mono_tints[3],
  #dropped_col,
  mono_tints[4],
  mono_tints[5]
  #dropped_col
)

mono_flow <- mon_plot_df %>%
  count(original, accepted, assessed, mapped, workflow)

mono <- flow_to_sankey(mono_flow)

mono_sankey <- sankeyNetwork(
  Links        = mono$links,
  Nodes        = mono$nodes,
  Source       = "source",
  Target       = "target",
  Value        = "value",
  NodeID       = "name",
  colourScale  = node_colour_js(mono_colours),
  nodeWidth    = 15,
  nodePadding  = 15,
  fontSize     = 13,
  fontFamily   = "Georgia",
  sinksRight   = TRUE,
  width        = 1000,
  height       = 300
)

# Make links grey by injecting a small JS tweak
# Re-apply the grey links JS, then save
mono_sankey <- htmlwidgets::onRender(mono_sankey, '
  function(el, x) {
    d3.select(el).selectAll(".link")
      .style("stroke", "#CCCCCC")
      .style("stroke-opacity", 0.5);

  }
')

mono_sankey

# ── SAVE AS PNG ──────────────────────────────────────────────────────────────
# Save as temp HTML then screenshot
saveWidget(mono_sankey, "02_figures/Figure_S2_sankey_monocots.html", selfcontained = TRUE)
webshot("02_figures/Figure_S2_sankey_monocots.html", "02_figures/Figure_S2_sankey_monocots.png", vwidth = 1000, vheight = 300, zoom = 2)

# ── LEGUMES ─────────────────────────────────────────────────────────────────
# Replace with your actual legume counts
leg_tints <- make_node_colours("#E69F00", 5)

leg_colours <- c(
  leg_tints[1],
  leg_tints[2],
  #dropped_col,
  leg_tints[3],
  #dropped_col,
  leg_tints[4],
  leg_tints[5]
  #dropped_col
)

leg_flow <- legs_plot_df %>%
  count(original, accepted, assessed, mapped, workflow)

legs <- flow_to_sankey(leg_flow)

legs_sankey <- sankeyNetwork(
  Links        = legs$links,
  Nodes        = legs$nodes,
  Source       = "source",
  Target       = "target",
  Value        = "value",
  NodeID       = "name",
  colourScale  = node_colour_js(leg_colours),
  nodeWidth    = 15,
  nodePadding  = 15,
  fontSize     = 13,
  fontFamily   = "Georgia",
  #sinksRight   = TRUE,
  width        = 1000,
  height       = 300
)

# Make links grey by injecting a small JS tweak
# Re-apply the grey links JS, then save
legs_sankey <- htmlwidgets::onRender(legs_sankey, '
  function(el, x) {
    d3.select(el).selectAll(".link")
      .style("stroke", "#CCCCCC")
      .style("stroke-opacity", 0.5);
  }
')

#d3.select(el).selectAll(".node text:last-child").remove();

legs_sankey

# ── SAVE AS PNG ──────────────────────────────────────────────────────────────
# Save as temp HTML then screenshot
saveWidget(legs_sankey, "02_figures/Figure_S2_sankey_legumes.html", selfcontained = TRUE)
webshot("02_figures/Figure_S2_sankey_legumes.html", "02_figures/Figure_S2_sankey_legumes.png", vwidth = 1000, vheight = 300, zoom = 2)


# # extra:
# # Load in the wcvp
# # Version 13 2024
# wcvp_names <- read_delim(paste0(work_dir,"00_raw_data/wcvp_names.csv"),
#                          delim = "|",
#                          escape_double = FALSE,
#                          trim_ws = TRUE)
#
# # get the number of accepted species
# wcvp_acc <-
#   wcvp_names |>
#   filter(taxon_rank == "Species",
#          taxon_status == "Accepted",
#          !is.na(taxon_name),
#          is.na(genus_hybrid),
#          is.na(species_hybrid)) |>
#   select(powo_id, taxon_name, family, genus, lifeform_description, climate_description)
# nrow(wcvp_acc)
#
# # get number of Monocots and Legumes species
# # Legumes - easy, just filter on Fabs
# legumes_all <- wcvp_acc %>%
#   dplyr::filter(family == "Fabaceae")
#
# legumes_all<- legumes_all %>% count()
#
# # monocots - need to get mon fams
# APGIV_families <- read.csv(paste0(work_dir,"00_raw_data/APGIV_families_orders_IUCN_crosswalk.csv"))
#
# monocots_all <- wcvp_acc %>%
#   left_join(APGIV_families, by = c(family = 'APG.IV.Family'))
#
# monocots_tot <- monocots_all %>% dplyr::filter(HigherGroups == "Monocots") %>% count()

