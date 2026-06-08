
# Non-genuine and genuine changes - tables and upset charts for comparison

library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggthemes)
library(ggupset)
library(gt)
library(webshot2)
library(cowplot)
library(magick)
library(rsvg)

# non - genuine changes ----
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load SRLI main table
SRLI_2026 <- read_csv(paste0(work_dir, "00_raw_data/SRLI_2026.csv")) %>%
  dplyr::filter(SRLI_2024_analysis == "include") %>%
  dplyr::filter(complete.cases(SRLI_2020_final_cat)) %>%
  dplyr::filter(SRLI_2020_final_cat != "NE")
#RLI_stat <- SRLI_2026 %>% filter(SRLI_2024_analysis == "include") # only the included species


# group and summarise changes 2000 - 2010
groups_2000_2010<- SRLI_2026 %>%
  group_by(group, change_2000_2010) %>% count()

# widen the data for the summary table
wide_2000_2010 <- groups_2000_2010 %>%
  pivot_wider(
    names_from = change_2000_2010,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(change_period = "change_2000_2010")

# group and summarise changes 2010 - 2020
groups_2010_2020<- SRLI_2026 %>%
  group_by(group, change_2010_2020) %>% count()

# widen the data for the summary table
wide_2010_2020 <- groups_2010_2020 %>%
  pivot_wider(
    names_from = change_2010_2020,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(change_period = "change_2010_2020")

# bind them together
change_summary <- bind_rows(wide_2000_2010,wide_2010_2020)

# add percentage
change_summary <- change_summary %>%
  rowwise() %>%
  mutate(
    total = sum(c_across(where(is.numeric)), na.rm = TRUE)
  ) %>%
  mutate(
    across(
      where(is.numeric) & !matches("^total$"),
      ~ round((.x / total) * 100,1),
      .names = "{.col}_pct"
    )
  ) %>%
  ungroup() %>%
  select(-total) %>%
  select(change_period,group,Genuine, Non_genuine, No_change, Genuine_pct, Non_genuine_pct, No_change_pct)


totals <- change_summary |>
  group_by(change_period) |>
  summarise(
    group       = "Total",
    Genuine     = sum(Genuine),
    Non_genuine = sum(Non_genuine),
    No_change   = sum(No_change),
    .groups = "drop"
  ) |>
  mutate(
    Genuine_pct     = Genuine     / (Genuine + No_change + Non_genuine) * 100,
    Non_genuine_pct = Non_genuine / (Genuine + No_change + Non_genuine) * 100,
    No_change_pct   = No_change   / (Genuine + No_change + Non_genuine) * 100

  )


# write this table to results and merge as a single excel file in sup mat
write.csv(change_summary,
          "02_figures/change_summary.csv",
          row.names = FALSE)

### formatting table for panel
library(dplyr)
library(gt)

my_table <- bind_rows(change_summary, totals) |>
  mutate(
    period = case_when(
      change_period == "change_2000_2010" ~ "2000–2010",
      change_period == "change_2010_2020" ~ "2010–2020"
    ),
    group = factor(group, levels = c("Legumes", "Monocots", "Total"))
  ) |>
  arrange(period, group) |>
  select(period, group, Genuine,Non_genuine,No_change,
         Genuine_pct,Non_genuine_pct,No_change_pct) |>
  group_by(period) |>
  gt(rowname_col = "group") |>

cols_add(
  Total = Genuine + Non_genuine + No_change
) |>

cols_label(
  Genuine          = "Genuine",
  Non_genuine      = "Non-genuine",
  No_change        = "No change",
  Total            = "Total",
  Genuine_pct      = "Genuine (%)",
  Non_genuine_pct  = "Non-genuine (%)",
  No_change_pct    = "No change (%)",
) |>

tab_spanner(
  label   = "Count",
  columns = c(Genuine, Non_genuine,No_change, Total)
) |>
  tab_spanner(
    label   = "Percentage",
    columns = c(Genuine_pct, Non_genuine_pct, No_change_pct)
  ) |>

fmt_integer(columns = c(Genuine, Non_genuine, No_change, Total)) |>
  fmt_number(
    columns  = c(Genuine_pct, Non_genuine_pct, No_change_pct),
    decimals = 1,
    pattern  = "{x}%"
  ) |>

tab_style(
    style     = cell_borders(sides = "top", color = "black", weight = px(1.5)),
    locations = cells_column_labels(
      columns = c(Genuine, Non_genuine, No_change, Total)
    )
  ) |>
    tab_style(
      style     = cell_borders(sides = "top", color = "black", weight = px(1.5)),
      locations = cells_column_labels(
        columns = c(Genuine_pct, Non_genuine_pct, No_change_pct)
      )
    ) |>
  tab_style(
    style     = cell_borders(sides = "all", color = "transparent"),
    locations = cells_stub(rows = group == "Total")
  ) |>
  tab_options(
    table.font.size                   = px(11),
    data_row.padding                  = px(0),
    row_group.padding                 = px(0),
    column_labels.padding             = px(0),
    column_labels.border.top.style    = "none",
    column_labels.border.bottom.style = "none",
    table_body.border.top.style       = "none",
    table.border.left.style           = "none",
    table.border.right.style          = "none",
    table_body.vlines.style           = "none",
    stub.border.style                 = "none"
  )

#gtsave(my_table, "table.pdf")
gtsave(my_table, "02_figures/Table_1.png", expand = 20, zoom = 3)

# add significance test? ----
# compare proportion of non genuine in 2000-2010 compared to 2010-2020
# same chi sq test?

# non - genuine - raw counts ----
# only for the non genuine (combined years, but split by group)
non_gen_2000_2010 <- SRLI_2026 %>% dplyr::filter(
  change_2000_2010 == "Non_genuine"
) %>%
  select(group, reason_2000_2010) %>%
  rename(combined_change_reason = reason_2000_2010)

non_gen_2010_2020 <- SRLI_2026 %>% dplyr::filter(
  change_2010_2020 == "Non_genuine"
) %>%
  select(group, reason_2010_2020) %>%
  rename(combined_change_reason = reason_2010_2020)

non_gen <- bind_rows(non_gen_2000_2010, non_gen_2010_2020)

# non-genuine changes reason
groups <- non_gen %>%
  group_by(combined_change_reason) %>%
  count()

write.csv(groups,
          "02_figures/Table_S5_non_genuine_raw_counts.csv",
          row.names = FALSE)

# chart
fit <- c("Knowledge of criteria" = 1,
               "New Information" = 116,
               "New Information&Incorrect data" = 5,
               "New Information&Incorrect data&Taxonomy" = 1,
               "New Information&Knowledge of criteria" = 4,
               "New Information&Knowledge of criteria&Criteria Revision" = 2,
               "New Information&Knowledge of criteria&Incorrect data" = 1,
               "New Information&Taxonomy" = 1,
               "Not reported" = 1)

# non - gen ggupset ----
non_gen_plot <- non_gen %>%
  mutate(
    # Split the semicolon-separated string into a list of vectors
    sets = strsplit(combined_change_reason, ";")
  )

up_non_gen <- ggplot(non_gen_plot, aes(x = sets, fill = group)) +
  geom_bar() +
  scale_x_upset() +
  scale_fill_manual(values = c("Legumes" = "#E69F00", "Monocots" = "#0072B2")
                    , name = NULL) +
  labs(
    title = "Reasons for non-genuine changes",
    x = NULL,
    y = "Count",
    fill = "Group"
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0.5
  ) +
  theme_base(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.80),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# save it
ggsave(
  "02_figures/non_gen_upset.png",
  width = 16,
  height = 9,
  units = "cm",
  plot = up_non_gen
)

# genuine counts ----
# chart ----
df <- read_excel("C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/00_raw_data/RLI_genuine_reasons_2.xlsx")

df_long <- df %>%
  pivot_longer(cols = c(Primary, Secondary),
               names_to = "type",
               values_to = "threat") %>%
  filter(!is.na(threat))

genuine_threat_reasons <- df_long %>%
  group_by(threat) %>%
  count() %>%
  arrange(desc(n))


write.csv(genuine_threat_reasons,
          "02_figures/Table_S4_genuine_threat_reasons.csv",
          row.names = FALSE)

# genuine upset ----
# 1. Collapse threats into a list column per species
upset_data_gen <- df_long %>%
  group_by(WCVP_accepted_name, group, powo_ID) %>%
  summarise(threats = list(unique(threat)), .groups = "drop")

up_gen <- ggplot(upset_data_gen, aes(x = threats, fill = group)) +
  geom_bar() +
  scale_x_upset() +
  scale_fill_manual(values = c("Legumes" = "#E69F00", "Monocots" = "#0072B2"),
                    name = NULL) +
  labs(
    title = "Reasons for genuine changes",
    x = NULL,
    y = "Count",
    fill = "Group"
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.panel.line.size = 0.5
  ) +
  theme_base(base_size = 12) +
  theme(
    legend.position = c(0.85, 0.80),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# save it
ggsave(
  "02_figures/gen_upset.png",
  width = 16,
  height = 9,
  units = "cm",
  plot = up_gen
)

# transitions matrices ----
# split by monocot and legumes
mono = dplyr::filter(SRLI_2026, group == "Monocots")
legs = dplyr::filter(SRLI_2026, group == "Legumes")

# Create custom alpha function - figure 3 matrix
alpha_transform <- function(x) {
  ifelse(is.na(x), NA,
         ifelse(x == 0, zero_alpha,
                min_alpha + (max_alpha - min_alpha) * (x - alpha_min) / (alpha_max - alpha_min)))
}

# monocots first
df <- data.frame(
  RL_2000 = mono$SRLI_2000_final_cat,
  RL_2020 = mono$SRLI_2020_final_cat)

status_levels <- c("LC", "NT", "VU", "EN", "CR") #"DD"
status_colors <- c(
  "LC" = "#60C659",
  "NT" = "#CCE226",
  "VU" = "#F9E814",
  "EN" = "#FC7F3F",
  "CR" = "#D81E05"
  #"DD" = "#D1D1C6"
)

# Count all transitions (using 2000 to 2020 as in your example)
counts <- df %>%
  count(From = RL_2000, To = RL_2020) %>%
  complete(From = status_levels, To = status_levels, fill = list(n = 0))

# set backgorun i.e. no change to grey
counts_v1 <- counts %>%
  mutate(
    FillColor = ifelse(From == To, "grey80", status_colors[To]),
    # Create alpha values that ensure count=1 is always visible
    AlphaVal = ifelse(From == To, NA, n)
  )

# Set manual scale parameters
alpha_min <- 1   # minimum meaningful count
alpha_max <- 5  # maximum count across all your datasets
zero_alpha <- 0.15  # alpha for zero counts (should be noticeably different from count=1)
min_alpha <- 0.4   # alpha for count=1 (should be clearly visible)
max_alpha <- 1.0   # alpha for maximum count

# plot it
mon_transitions <- ggplot(counts_v1, aes(x = From, y = To, fill = FillColor)) +
  geom_tile(aes(alpha = alpha_transform(AlphaVal)), color = "white", size = 1) +
  geom_text(aes(label = n), size = 5, color = "black", fontface = "bold") +
  scale_fill_identity() +
  scale_alpha_identity() +  # Use identity since we're manually transforming
  scale_x_discrete(limits = status_levels, position = "top") +
  scale_y_discrete(limits = rev(status_levels)) +
  theme_base(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(
    x = "2000 Red List category",
    y = "2020 Red List category"
  )

ggsave("02_figures/monocots2000-2020_grid.png", width = 6, height = 6, plot = mon_transitions)

# legumes next
df <- data.frame(
  RL_2000 = legs$SRLI_2000_final_cat,
  RL_2020 = legs$SRLI_2020_final_cat)

# status levels and colours as for monocots

# Count all transitions (using 2000 to 2020 as in your example)
counts <- df %>%
  count(From = RL_2000, To = RL_2020) %>%
  complete(From = status_levels, To = status_levels, fill = list(n = 0))

# set backgorun i.e. no change to grey
counts_v1 <- counts %>%
  mutate(
    FillColor = ifelse(From == To, "grey80", status_colors[To]),
    # Create alpha values that ensure count=1 is always visible
    AlphaVal = ifelse(From == To, NA, n)
  )

# Set manual scale parameters
alpha_min <- 1   # minimum meaningful count
alpha_max <- 5  # maximum count across all your datasets
zero_alpha <- 0.15  # alpha for zero counts (should be noticeably different from count=1)
min_alpha <- 0.4   # alpha for count=1 (should be clearly visible)
max_alpha <- 1.0   # alpha for maximum count

# plot it
leg_transitions <- ggplot(counts_v1, aes(x = From, y = To, fill = FillColor)) +
  geom_tile(aes(alpha = alpha_transform(AlphaVal)), color = "white", size = 1) +
  geom_text(aes(label = n), size = 5, color = "black", fontface = "bold") +
  scale_fill_identity() +
  scale_alpha_identity() +  # Use identity since we're manually transforming
  scale_x_discrete(limits = status_levels, position = "top") +
  scale_y_discrete(limits = rev(status_levels)) +
  theme_base(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(
    x = "2000 Red List category",
    y = "2020 Red List category"
  )

ggsave("02_figures/legumes2000-2020_grid.png", width = 6, height = 6,
       plot = leg_transitions)

# combined panel ----

# Combine with labels
add_border <- function(plot, colour) {
  ggdraw(plot) +
    theme(plot.background = element_rect(colour = colour, fill = NA, linewidth = 2))
}

add_padding <- function(plot, pad = 0.05) {
  ggdraw() +
    draw_plot(plot, x = pad, y = pad, width = 1 - 2*pad, height = 1 - 2*pad)
}

recolour_svg <- function(svg_path, colour, size = 200) {
  image_read_svg(svg_path, width = size) |>
    image_colorize(opacity = 100, color = colour)
}

pad <- 0.1

mon_transitions_bordered <- add_border(mon_transitions, "#0072B2")
leg_transitions_bordered <- add_border(leg_transitions, "#E69F00")

leg_icon <- recolour_svg("02_figures/PhyloPic_Legume.svg", "#E69F00")
mon_icon <- recolour_svg("02_figures/PhyloPic_Anomochloa_marantoidea.svg", "#0072B2")

mon_transitions_final <- ggdraw() +
  draw_plot(mon_transitions_bordered, x = pad, y = pad, width = 1 - pad, height = 1 - 2*pad) +
  draw_image(mon_icon, x = 0.01, y = 0.82, width = 0.15, height = 0.15)

leg_transitions_final <- ggdraw() +
  draw_plot(leg_transitions_bordered, x = pad, y = pad, width = 1 - pad, height = 1 - 2*pad) +
  draw_image(leg_icon, x = 0.01, y = 0.82, width = 0.15, height = 0.15)

up_gen_final <- up_gen + theme(plot.margin = margin(l = 60, t = 5, r = 5, b = 5))
up_non_gen_final <- up_non_gen + theme(plot.margin = margin(l = 60, t = 5, r = 5, b = 5))

top_row <- plot_grid(
  up_gen, up_non_gen, mon_transitions_final, leg_transitions_final,
  labels = c("A", "B", "C", "D"),
  nrow   = 2
)
top_row

ggsave("02_figures/Figure_2_ABCD.png", top_row, width = 13, height = 10)




# Extras
# fit <- euler(c("Knowledge of criteria" = 1,
#                "New Information" = 116,
#                "New Information&Incorrect data" = 5,
#                "New Information&Incorrect data&Taxonomy" = 1,
#                "New Information&Knowledge of criteria" = 4,
#                "New Information&Knowledge of criteria&Criteria Revision" = 2,
#                "New Information&Knowledge of criteria&Incorrect data" = 1,
#                "New Information&Taxonomy" = 1,
#                "Not reported" = 1),
#              shape = "ellipse")
#
# error_plot(fit)
#
# png(
#   file   = "02_figures/non_gen_Euler.png",
#   width  = 10,      # pixels
#   height = 5,     # pixels
#   res    = 300,      # DPI (resolution)
#   units  = "in"      # or "in" if you want width/height in inches
# )
#
# plot(fit,
#      #quantities = list(type = c("percent")),
#      fills = list(fill = c("cornsilk",
#                            "lightblue1",
#                            "darkorange",
#                            "darkgrey",
#                            "plum2")),
#      legend = list(side = "right")
# )
#
# dev.off()

# # Create an UpsetR Plot
# png(
#   file   = "02_figures/non_gen_upset.png",
#   width  = 6,      # pixels
#   height = 5,     # pixels
#   res    = 600,      # DPI (resolution)
#   units  = "in"      # or "in" if you want width/height in inches
# )
#
# upset(fromExpression(fit), order.by = "freq")
# queries = list(
#   list(
#     query = elements,
#     params = list("physicalActivityPerMonth", 0,4),
#     color = "#Df5286",
#     active = T,
#     query.name = "Physically Active 1x/Week or Less"
#   )
# )
# grid.text(
#   x = 0.90,
#   y = 0.05,
#   gp = gpar(
#     fontsize = 18,
#     fontface = 6
#   )
# )
# dev.off()

