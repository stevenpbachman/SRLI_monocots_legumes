
# Figure 4 ----

# Validation results
# stats test - gbif vs curated
# A - distance from aoh cell to validation point
# B - accuracy
# C - point prevalance
# stats - summary total species validated by group and method

# Comparison of TRUE/FALSE for each buffer. Bar chart

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
library(cowplot)

# organise data
# validation stat test 1 ----
# gbif vs curated
# start from scratch so that this is standalone
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"
SRLI_2024 <- read_csv(paste0(work_dir,"00_raw_data/SRLI_2024.csv"))
RLI_stat <- SRLI_2024 %>% filter(SRLI_2024_analysis == "include") # only the included species
RLI_stat_thin <- RLI_stat %>% select(group, WCVP_current_id) # only ID and group

# read in the validation summary files
val_gbif <- list.files(path=paste0(work_dir,"/04_outputs/gbif_workflow/validation_summary/"),
                       pattern="csv",
                       all.files=FALSE,
                       full.names=TRUE,
                       recursive=TRUE)

gbif_read_in <- map_df(val_gbif, read_csv, show_col_types = FALSE)
gbif_read_in$source <- "GBIF"

val_curated <- list.files(path=paste0(work_dir,"/04_outputs/curated_workflow/validation_summary/"),
                          pattern="csv",
                          all.files=FALSE,
                          full.names=TRUE,
                          recursive=TRUE)

curated_read_in <- map_df(val_curated, read_csv, show_col_types = FALSE)
curated_read_in$source <- "Curated"

# combine the gbif and curated results
comb_val <- dplyr::bind_rows(curated_read_in, gbif_read_in)
comb_val <- comb_val %>% left_join(RLI_stat_thin, by = c("sp_id" = "WCVP_current_id"))

# set factors to get order and friendly labels
comb_val$group <- factor(comb_val$group, levels = c("Monocots", "Legumes"))
comb_val$method <- str_remove(comb_val$method, "_buffer")
comb_val$method <- factor(comb_val$method, levels = c("exact", "1km", "2km", "5km"))
comb_val$source <- str_replace(comb_val$source, "gbif", "GBIF")
#comb_val$source <- str_replace(comb_val$source, "curated", "Curated")

comb_val <- comb_val %>% filter(complete.cases(method))
#comb_val$method <- as.character(comb_val$method)
exact <- comb_val %>% dplyr::filter(method == "exact")
buff_1km <- comb_val %>% dplyr::filter(method == "1km")
buff_2km <- comb_val %>% dplyr::filter(method == "2km")
buff_5km <- comb_val %>% dplyr::filter(method == "5km")

# Stats test - Create a contingency table - exact
contingency_table_exact <- table(exact$source, exact$valid_test)
contingency_table_buff_1km <- table(buff_1km$source, buff_1km$valid_test)
contingency_table_buff_2km <- table(buff_2km$source, buff_2km$valid_test)
contingency_table_buff_5km <- table(buff_5km$source, buff_5km$valid_test)

# Perform Fishers exact Test
fisher_test_exact <- fisher.test(contingency_table_exact)
print(fisher_test_exact)

fisher_test_buff_1km <- fisher.test(contingency_table_buff_1km)
print(fisher_test_buff_1km)

fisher_test_buff_2km <- fisher.test(contingency_table_buff_2km)
print(fisher_test_buff_2km)

fisher_test_buff_5km <- fisher.test(contingency_table_buff_5km)
print(fisher_test_buff_5km)

# explain the results:
# exact - curated is better, but by 2km there is hardly any difference, not significant
#invert the odds ratio: e.g  1 / 0.38 = 2.63

# 4A ----
#summarize first: calculate % of species better than random
val_summary <- comb_val %>%
  group_by(source, group, method) %>%
  summarise(
    success_rate = mean(valid_test == TRUE),
    .groups = "drop"
  )

# Plot proportions directly
buffer_comparison <- ggplot(val_summary, aes(x = method, y = success_rate, fill = source)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "grey30") +
  facet_wrap(~ group) +
  geom_hline(yintercept = 0.95, color = "black", linetype = "dashed", size = 0.5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("Curated" = "#CC79A7", "GBIF" = "#509E2F")) +
  labs(
    x = "Buffer radius from points",
    y = "Proportion better than random",
    fill = "Group"
  ) +
  theme_base(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ),
  strip.text = element_text(
    size = 14,
    face = "bold"
  ),
    #strip.background = element_rect(fill = NA, colour = "grey", size = 0.5), # lighter, thinner outline
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

buffer_comparison
ggsave("02_figures/Figure_4_A.png", width = 10, height = 4)

# 4B ----
# now get the accuracy values from comb_val
comb_val$diff <- comb_val$meanPointsprev - comb_val$meanAOHprev

acc <- comb_val %>%
  group_by(group, method, source) %>%
  summarise(
    #n = n(),
    median_accuracy = median(diff),
    iqr_accuracy = IQR(diff),
    .groups = "drop"
  )

acc_box <- ggplot(comb_val, aes(x = method, y = diff, fill = source)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.9) +
  scale_fill_manual(values = c("Curated" = "#CC79A7", "GBIF" = "#509E2F")) +
  facet_wrap(~ group) +
  labs(
    x = "Validation method",
    y = "Accuracy",
    fill = "Source"
  ) +
  theme_base(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = 'NA'),
    strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ),
  strip.text = element_text(
    size = 14,
    face = "bold"
  ),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

acc_box
ggsave("02_figures/Figure_6_B.png", width = 10, height = 4)

# 4C ----
#prevalence validation plot
# monocots and legumes combined left and right
combo_1km <- comb_val %>% dplyr::filter(method == "1km")
combo_1km$diff <- combo_1km$meanPointsprev - combo_1km$meanAOHprev

# plot
combo_1km_plot <- ggplot(
  data = combo_1km,
  aes(x = meanAOHprev, y = meanPointsprev)
) +
  geom_point(
    aes(
      size = sumpointstot,
      fill = source,        # now color points by source instead of group
      alpha = abs(diff)
    ),
    shape = 21,
    color = "grey22",
    stroke = 0.2
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey22", linewidth = 0.25) +
  scale_fill_manual(
    name = "Source",
    values = c("Curated" = "#CC79A7", "GBIF" = "#509E2F")
  ) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  scale_size_continuous(name = "Number of points", range = c(1, 8)) +
  facet_grid(~ group) +  # still facet by group
  labs(
    x = "Mean AOH prevalence",
    y = "Mean Point Prevalence"
  ) +
  scale_x_continuous(
    breaks = c(0, 1),
    labels = c("0", "1"),
    limits = c(0, 1)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(
    alpha = "none",
    fill = guide_legend(override.aes = list(size = 5))
  ) +
  theme_base(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = 'NA'),
    strip.background = element_rect(
    fill = "white",
    colour = "black",
    linewidth = 0.5
  ),
  strip.text = element_text(
    size = 14,
    face = "bold"
  ),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

combo_1km_plot

ggsave("02_figures/Figure_4_B.png", width = 10, height = 4)

# combine in single figure
# Figure 6 ----
combo <- cowplot::plot_grid(buffer_comparison, acc_box, combo_1km_plot, nrow = 3, labels = c('A', 'B', 'C'))
ggsave("02_figures/Figure_4_validation_results_ABC.png", width = 10, height = 12)

# validation stat test 2 ----
# summary total species validated by group and method

# comb_val includes all validation results
aoh_val_tested <- comb_val %>% group_by(method) %>% count()
print(aoh_val_tested$n[1])

# read in errors and summarise
# GBIF
aoh_errors_gbif <- list.files(path=paste0(work_dir,"04_outputs/gbif_workflow/errors/"),
                              pattern="csv",
                              all.files=FALSE,
                              full.names=TRUE,
                              recursive=TRUE)

# read in the AOH error files
gbif_errors_read_in = purrr::map_df(aoh_errors_gbif, read_csv, show_col_types = FALSE)

# now get the errors
# curated
aoh_errors_curated <- list.files(path=paste0(work_dir,"04_outputs/curated_workflow/errors/"),
                                 pattern="csv",
                                 all.files=FALSE,
                                 full.names=TRUE,
                                 recursive=TRUE)

# read in the AOH error files
curated_aoh_read_in = purrr::map_df(aoh_errors_curated, read_csv, show_col_types = FALSE) %>%
  dplyr::mutate(source = "curated")

combined_errors <- bind_rows(curated_aoh_read_in, gbif_errors_read_in)

sum_error_reasons <- combined_errors %>% group_by(error) %>% count()

# Table S7 ----
# make a table with the proportions - sup mat
comb_val$diff <- comb_val$meanPointsprev - comb_val$meanAOHprev

grouped_summary <- comb_val %>%
  group_by(group, method, source) %>%
  summarise(proportion_TRUE = mean(valid_test == TRUE),
            mean_accuracy = mean(diff, na.rm = TRUE),
            sd_accuracy = sd(diff, na.rm = TRUE),
            count = n())

# Summary without group (just method and source)
method_source_summary <- comb_val %>%
  group_by(method, source) %>%
  summarise(proportion_TRUE = mean(valid_test == TRUE),
            mean_accuracy = mean(diff, na.rm = TRUE),
            sd_accuracy = sd(diff, na.rm = TRUE),
            count = n(),
            group = "All") # Add a placeholder for group

# Summary without group (just method )
method_summary <- comb_val %>%
  group_by(method) %>%
  summarise(proportion_TRUE = mean(valid_test == TRUE),
            mean_accuracy = mean(diff, na.rm = TRUE),
            sd_accuracy = sd(diff, na.rm = TRUE),
            count = n(),
            group = "All") # Add a placeholder for group

# Combine them
Table_S7_validation_proportions <- bind_rows(grouped_summary,
                                             method_source_summary,
                                             method_summary)


write.csv(Table_S7_validation_proportions,
          "02_figures/Table_S7_validation_proportions.csv",
          row.names = FALSE)


# supplementary validation figures

# Figure S5 validation combo ----
# facet wrap version to compare different buffers

# 1 GBIF legumes
comb_val$diff <- comb_val$meanPointsprev - comb_val$meanAOHprev

gbif_val_legs <- comb_val %>% dplyr::filter(source == "GBIF", group == "Legumes")

gbif_val_legs_comparison <- ggplot(data = gbif_val_legs, aes(x = meanAOHprev, y = meanPointsprev)) +
  geom_point(aes(size = sumpointstot, fill = diff), #diff
             shape = 21, color = "grey22", stroke = 0.2) +  # Black outline and fill based on 'diff'
  geom_abline(slope = 1, intercept = 0, color = "grey22", linewidth = 0.25) +
  scale_fill_gradient2(low = "white", high = "#E69F00",
                       midpoint = 0, name = "Points - AOH difference") +  # Color gradient for 'diff'
  scale_size_continuous(name = "Number of points", range = c(1, 8)) +  # Control point sizes
  theme_base(base_size = 16) +
  xlim(0, 1) + ylim(0, 1) +
  facet_grid(~ group ~ method) +
  labs(x = "Mean AOH prevalence", y = "Mean Point Prevalence") +
  theme_base(base_size = 14) +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  scale_x_continuous(
    breaks = c(0, 1),  # Only show ticks at 0 and 1
    labels = c("0", "1"),  # Corresponding labels for those breaks
    limits = c(0, 1)  # Set the axis limits
  ) +
  scale_y_continuous(limits = c(0, 1))

gbif_val_legs_comparison

# 2 GBIF monocots
gbif_val_mons <- comb_val %>% dplyr::filter(source == "GBIF", group == "Monocots")

gbif_val_mons_comparison <- ggplot(data = gbif_val_mons, aes(x = meanAOHprev, y = meanPointsprev)) +
  geom_point(aes(size = sumpointstot, fill = diff), #diff
             shape = 21, color = "grey22", stroke = 0.25) +  # Black outline and fill based on 'diff'
  geom_abline(slope = 1, intercept = 0, color = "grey22", linewidth = 0.1) +
  scale_fill_gradient2(low = "white", high = "#0072B2",
                       midpoint = 0, name = "Points - AOH difference") +  # Color gradient for 'diff'
  scale_size_continuous(name = "Number of points", range = c(1, 8)) +  # Control point sizes
  theme_base(base_size = 16) +
  #xlim(0, 1) + ylim(0, 1) +
  facet_grid(~ group ~ method) +
  labs(x = "Mean AOH prevalence", y = "Mean Point Prevalence") +
  theme_base(base_size = 14) +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  # Customize x-axis to only show ticks at 0 and 1
  scale_x_continuous(
    breaks = c(0, 1),  # Only show ticks at 0 and 1
    labels = c("0", "1"),  # Corresponding labels for those breaks
    limits = c(0, 1)  # Set the axis limits
  ) +
  scale_y_continuous(limits = c(0, 1))

gbif_val_mons_comparison

# 3 Curated legumes
curated_val_legs <- comb_val %>% dplyr::filter(source == "Curated", group == "Legumes")

curated_val_legs_comparison <- ggplot(data = curated_val_legs, aes(x = meanAOHprev, y = meanPointsprev)) +
  geom_point(aes(size = sumpointstot, fill = diff), #diff
             shape = 21, color = "grey22", stroke = 0.2) +  # Black outline and fill based on 'diff'
  geom_abline(slope = 1, intercept = 0, color = "grey22", linewidth = 0.25) +
  scale_fill_gradient2(low = "white", high = "#E69F00",
                       midpoint = 0, name = "Points - AOH difference") +  # Color gradient for 'diff'
  scale_size_continuous(name = "Number of points", range = c(1, 8)) +  # Control point sizes
  theme_base(base_size = 16) +
  xlim(0, 1) + ylim(0, 1) +
  facet_grid(~ group ~ method) +
  labs(x = "Mean AOH prevalence", y = "Mean Point Prevalence") +
  theme_base(base_size = 14) +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  scale_x_continuous(
    breaks = c(0, 1),  # Only show ticks at 0 and 1
    labels = c("0", "1"),  # Corresponding labels for those breaks
    limits = c(0, 1)  # Set the axis limits
  ) +
  scale_y_continuous(limits = c(0, 1))

curated_val_legs_comparison

# 4 Curated monocots
curated_val_mons <- comb_val %>% dplyr::filter(source == "Curated", group == "Monocots")

curated_val_mons_comparison <- ggplot(data = curated_val_mons, aes(x = meanAOHprev, y = meanPointsprev)) +
  geom_point(aes(size = sumpointstot, fill = diff), #diff
             shape = 21, color = "grey22", stroke = 0.2) +  # Black outline and fill based on 'diff'
  geom_abline(slope = 1, intercept = 0, color = "grey22", linewidth = 0.25) +
  scale_fill_gradient2(low = "white", high = "#0072B2",
                       midpoint = 0, name = "Points - AOH difference") +  # Color gradient for 'diff'
  scale_size_continuous(name = "Number of points", range = c(1, 8)) +  # Control point sizes
  theme_base(base_size = 16) +
  xlim(0, 1) + ylim(0, 1) +
  facet_grid(~ group ~ method) +
  labs(x = "Mean AOH prevalence", y = "Mean Point Prevalence") +
  theme_base(base_size = 14) +
  theme(
    legend.text = element_text(size = 10),
    legend.position = "none",
    legend.background = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  scale_x_continuous(
    breaks = c(0, 1),  # Only show ticks at 0 and 1
    labels = c("0", "1"),  # Corresponding labels for those breaks
    limits = c(0, 1)  # Set the axis limits
  ) +
  scale_y_continuous(limits = c(0, 1))

curated_val_mons_comparison

# combination plot
combo <- plot_grid(gbif_val_mons_comparison, gbif_val_legs_comparison,
                   curated_val_mons_comparison, curated_val_legs_comparison,
                   nrow = 4, labels = c('A', 'B', 'C', 'D'))

ggsave("02_figures/Figure_S4_validation_combo.png", width = 14, height = 15)


