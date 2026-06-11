
# Red list index trend lines ----

# Figure 3 - Red List Index

# set up working directory
#work_dir = "G:/Users/sb42kg/SRLI_2024/"
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load the function
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

# load relevant files
SRLI_2024 <- read_csv(paste0(work_dir,"00_raw_data/SRLI_2026.csv"))
RLI_stat <- SRLI_2024 %>% filter(SRLI_2024_analysis == "include") # only the included species

# split by monocot and legumes
mono = dplyr::filter(RLI_stat, group == "Monocots")
nrow(mono)
legs = dplyr::filter(RLI_stat, group == "Legumes")
nrow(legs)

# summarise data to get ready for ggplot
rli_mono_2000 <- table(mono$SRLI_2000_final_cat)
rli_mono_2010 <- table(mono$SRLI_2010_final_cat)
rli_mono_2020 <- table(mono$SRLI_2020_final_cat)
rli_legs_2000 <- table(legs$SRLI_2000_final_cat)
rli_legs_2010 <- table(legs$SRLI_2010_final_cat)
rli_legs_2020 <- table(legs$SRLI_2020_final_cat)

# all RLI values - bootstrap
# final categories combined
rli_boot_m2000 <- rli_boot(rli_mono_2000, boot = TRUE, year = 2000, group = "Monocots", dd_species = unname(rli_mono_2000["DD"]))
rli_boot_m2010 <- rli_boot(rli_mono_2010, boot = TRUE, year = 2010, group = "Monocots", dd_species = unname(rli_mono_2010["DD"]))
rli_boot_m2020 <- rli_boot(rli_mono_2020, boot = TRUE, year = 2020, group = "Monocots", dd_species = unname(rli_mono_2020["DD"]))

rli_boot_l2000 <- rli_boot(rli_legs_2000, boot = TRUE, year = 2000, group = "Legumes", dd_species = unname(rli_legs_2000["DD"]))
rli_boot_l2010 <- rli_boot(rli_legs_2010, boot = TRUE, year = 2010, group = "Legumes", dd_species = unname(rli_legs_2010["DD"]))
rli_boot_l2020 <- rli_boot(rli_legs_2020, boot = TRUE, year = 2020, group = "Legumes", dd_species = unname(rli_legs_2020["DD"]))

# join them up
rli_all <- bind_rows(rli_boot_m2000, rli_boot_m2010, rli_boot_m2020,
                     rli_boot_l2000, rli_boot_l2010, rli_boot_l2020)

# Table S3 - RLI values tabular----
# save down the table so we can add to sup mat
Table_S6_rli_values <- rli_all

write.csv(Table_S6_rli_values,
          "02_figures/Table_S6_rli_values.csv",
          row.names = FALSE)

#Result 5. ----
Table_S6_rli_values %>%
  filter(Year %in% c("2000", "2020")) %>%
  arrange(Group, Year) %>%
  group_by(Group) %>%
  summarise(text = glue("Result 5
                        {first(Group)} RLI went from {round(first(RLI), 3)}% in 2000 to {round(last(RLI), 3)}% in 2020")) %>%
  pull(text)


# calculate equivalent change in categories based on RLI value diffs
# Subtract the 2000 RLI value from the 2020 value – this is the RLI change value
# A decline of 1/5 = 0.20 in the RLI would correspond to every species moving up
# exactly one full category (e.g., all LC → NT, or all VU → EN).
# so divide the drop by 0.2
# That should give you the percentage of species in your dataset (non DD) moving
# one red-list category closer to extinction.

# calculation
result <- Table_S6_rli_values %>%
  filter(Year %in% c(2000, 2020)) %>%              # only 2000 and 2020
  select(Group, Year, RLI) %>%
  pivot_wider(names_from = Year, values_from = RLI) %>%  # wide format
  mutate(
    diff = abs(`2020` - `2000`),                        # difference
    frac = diff / 0.2   # divide by 0.2
  )

result %>%
  mutate(
    perc = paste0(round(diff, 2), "%")       # keep 2 decimal
  )

#Result 6. ----
result %>%
  glue_data("Result 6
            {Group} fraction of species changing one category: {round((frac * 100),2)}%")

# for the h line - baseline value
y_int_m <- as.numeric(rli_all %>% filter(Year == 2000,
                                         Group == "Monocots") %>%
                        select(RLI))

y_int_l <- as.numeric(rli_all %>% filter(Year == 2000,
                                         Group == "Legumes") %>%
                        select(RLI))

# ggplot - SRLI only
RLI_simple = ggplot(rli_all, aes(x = Year, y = RLI, group = Group)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = Group), fill = "grey", alpha = 0.2) +
  #geom_hline(yintercept = c(0.8621901, 0.9019895),
  geom_hline(yintercept = c(y_int_m, y_int_l),
             color = "grey10",    # dark grey, near black
             linetype = "dotted",
             size = 0.8) +
  geom_line(aes(y = RLI, linetype = Group, color = Group), size = 1, alpha = 0.9) +
  geom_point(aes(y = RLI, color = Group), size = 3, alpha = 0.9) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00")) +
  scale_linetype_manual(values = c("Monocots" = "solid", "Legumes" = "solid")) +
  labs(
    y = "Red List Index of Species Survival",
    color = NULL,
    linetype = NULL
  ) +
  theme_base(base_size = 14) +
  theme(legend.position="bottom",
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA)) +
  ylim(0.85, 0.92) +
  scale_x_continuous(n.breaks=3)

# remove legend
RLI_simple <- RLI_simple +
  theme(legend.position = "none")

# save it
#ggsave("02_figures/Figure_3_RLI_simple.png", width = 4.5, height = 6)

# Make this part B of figure 2.
# now add in birds, mammals aggregated
RLI_others <- read_csv(paste0(work_dir,"00_raw_data/global_aggregated_RLI.csv")) #data to make the RLI plots for other groups

RLI_comb <- bind_rows(RLI_others, rli_all)

# the plot
RLI_all = ggplot(RLI_comb, aes(x = Year, y = RLI , group = Group, ymin=ci_lower, ymax=ci_upper)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, group = Group), fill = "grey", alpha = 0.2) +
  geom_line(aes(y = RLI, linetype = Group, color = Group), size = 1) +
  scale_color_manual(values = c("Monocots" = "#0072B2", "Legumes" = "#E69F00",
                                "Global aggregated" = "grey")) +
  scale_linetype_manual(values = c("Monocots" = "solid", "Legumes" = "solid",
                                   "Global aggregated" = "solid")) +
  geom_point(
    data = subset(RLI_comb, Group != "Global aggregated"),
    aes(y = RLI, color = Group),
    size = 3
  ) +
  labs(
    y = "Red List Index of Species Survival",
    color = NULL,
    linetype = NULL
  ) +
  theme_base(base_size = 14) +
  theme(legend.position="bottom",
        plot.background = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA)) +
  ylim(0.69, 0.92) +
  scale_x_continuous(n.breaks=3)

# remove legend and y axis
RLI_all_clean <- RLI_all +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = "none")

# joint figure
library(patchwork)
RLI_simple_clean <- RLI_simple +
  guides(colour = "none", fill = "none", linetype = "none")

combined <- RLI_simple_clean + RLI_all_clean +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.justification = "center"
  )

#leg_icon <- recolour_svg("02_figures/PhyloPic_Legume.svg", "#E69F00")
#mon_icon <- recolour_svg("02_figures/PhyloPic_Anomochloa_marantoidea.svg", "#0072B2")

# save it
ggsave("02_figures/Figure_3_RLI_aggregate.png", width = 7, height = 5.5)


# extras
# species numbers for each group
# Legumes <- RL_categories_per_group %>%
#   filter(group == "Legumes",
#          dataset == "SRLI 2020",
#          year == 2020,
#          category != "DD",
#          category != "NE") %>%
#   summarise(total = sum(count, na.rm = TRUE))
#
# Monocots <- RL_categories_per_group %>%
#   filter(group == "Monocots",
#          dataset == "SRLI 2020",
#          year == 2020,
#          category != "DD",
#          category != "NE") %>%
#   summarise(total = sum(count, na.rm = TRUE))
#
# species_counts <- data.frame("Legumes" = Legumes$total,
#                              "Monocots" = Monocots$total)
