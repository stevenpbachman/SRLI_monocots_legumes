
# Red List Index shows increasing extinction risk for Monocot and Legume plant species

## Overview

This repository contains R code used to perform the analysis described in the manuscript "Red List Index shows increasing extinction risk for Monocot and Legume plant species". 

## AOH map generation

The first part contains code used to generate area of habitat (AOH) maps. Note that the maps were analysed to generate results, but individual maps were not saved. 

The `analysis/XXX.R` script can be used to replicate the AOH workflow for curated data and `analysis/YYY.R` can be used to replicate the AOH workflow for the GBIF data.

## Analysis and figure generation

The second part outlines scripts used to analyse the data and generate the figures.
Start with `analysis/00_make_figures.R`

## Software 

reference the packages used

``` r
# install.packages("devtools")
devtools::install_github("stevenpbachman/LCr")
```

## Usage

Primary functions:

1.  `get_name_keys` check name status against taxonomic name backbones
2.  `get_gbif_occs` gather occurrence records from GBIF 
3.  `check_occs` & `clean_occs`  run occurrence quality checks and clean
4.  `make_sis_occs` generate the SIS point files
5.  `make_sis_csvs` generate SIS connect csv files
6.  `make_zip` make a zip file for SIS Connect imports

Note the following example is for plant species. For help on processing species in other kingdoms see the vignette [Fungi workflow](https://stevenpbachman.github.io/LCr/articles/SIS_connect_fungi_workflow.html).

### Example

``` r
library(LCr)

# create a species list
sp <- data.frame(sp_col = c("Crabbea acaulis", "Crabbea nana"))

# check names against GBIF and WCVP name backbones 
sp_keys <- LCr::get_name_keys(df = sp, name_column = "sp_col")

# get the raw GBIF occs and native ranges from Plants of the World Online
sp_occs <-  LCr::get_gbif_occs(sp_keys)
native_ranges <- LCr::get_native_range(keys = sp_keys)

# run the occurrence quality checks and clean data
checked_occs <- LCr::check_occs(sp_occs$points, native_ranges, buffer = 2000) # check
flagged_occs <- checked_occs$checked_data # flag
flagged_sum <- checked_occs$summary # summarise
cleaned_result <- LCr::clean_occs(flagged_occs) # clean

# define parameters for output files i.e. the person who is assessing the species
first_name <- "Steven"
second_name <- "Bachman"
email = "s.bachman@kew.org"
institution = "Royal Botanic Gardens, Kew"

# generate the point file
sis_point_file <- LCr::make_sis_occs(cleaned_result$clean_occs,
                                     first_name = first_name,
                                     second_name = second_name,
                                     institution = institution)

# generate output csv files
lc_sis_files <- LCr::make_sis_csvs(unique_id = sp_keys$GBIF_usageKey,
                                   wcvp_ipni_id = sp_keys$wcvp_ipni_id,
                                   first_name = "Steven",
                                   second_name = "Bachman",
                                   email = "s.bachman@kew.org",
                                   institution = "Royal Botanic Gardens, Kew",
                                   family = sp_keys$GBIF_family,
                                   genus = sp_keys$GBIF_genus,
                                   species = sp_keys$GBIF_species,
                                   gbif_ref = sp_occs,
                                   powo_ref = TRUE, 
                                   taxonomicAuthority = sp_keys$GBIF_authorship,
                                   occs = sis_point_file,
                                   kingdom = "plantae",
                                   native_ranges = native_ranges
)

# make the ZIP file to send to IUCN SIS Connect
make_zip(lc_sis_files)
```
### Citation

To cite `LCr` please run:

    citation("LCr")


