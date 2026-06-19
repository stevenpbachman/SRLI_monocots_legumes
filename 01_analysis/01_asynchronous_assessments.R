
# Asynchronous assessments ----

# show asynchronous RL assessments for SRLI species
# illustrates drift from assessment year 2010 due to numerous reasons
# assessments picked up by other projects e.g. 2017,2018 - global tree assessment
# makes it hard to show change over time between two time points, hence need
# to realign to reference years e.g. 2000,2010,2020

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)

work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load relevant files
SRLI_2026 <- read_csv(paste0(work_dir, "00_raw_data/SRLI_2026.csv")) %>%
  dplyr::filter(SRLI_2024_analysis == "include") %>%
  dplyr::filter(complete.cases(sis_assessments))

# organise
SRLI_2026 <- SRLI_2026 %>%
  separate(
    sis_assessments,
    into = c("year_start", "year_end"),
    sep = "-",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(across(c(year_start, year_end), as.integer))

# summarise
year_counts <- SRLI_2026 %>%
  pivot_longer(
    cols = c(year_start, year_end),
    names_to = "type",
    values_to = "year"
  ) %>%
  filter(!is.na(year), year >= 2005, year <= 2025) %>%
  count(year, group)

# plot it
ggplot(year_counts, aes(x = year, y = n, fill = group)) +
  geom_col() +
  scale_fill_manual(values = c("Monocots" = "#0072B2", "Legumes"  = "#E69F00"), name = NULL) +
  geom_vline(
    xintercept = 2010,
    linetype = "dotted",
    color = "black",
    size = 1
  ) +
  geom_vline(
    xintercept = 2020,
    linetype = "dotted",
    color = "black",
    size = 1
  ) +
  scale_x_continuous(labels = as.character(year_counts$year), breaks = year_counts$year) +
  ylab("Number of assessments") +
  xlab("Assessment year") +
  theme_base(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.85),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# save it
ggsave(
  "02_figures/Figure_S1_SRLI_assessment_years.png",
  width = 30,
  height = 14,
  units = "cm"
)
