

#' Utility functions needed for the analysis - source this first
#'
#' functions include:
#' AOH single - run a single AOH calcualtion for a stack e.g. 1992 - 2020
#' AOH stats - generate statistics based on AOH e.g. EOO and AOO from cells

# libraries ----
library(terra)
library(dplyr)
library(sf)
library(readr)
library(rCAT)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

# set up working directory
#work_dir = "G:/Users/sb42kg/SRLI_2024/"
work_dir = "C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/01_Projects/SRLI_project/SRLI_2_paper/SRLI_2024_C_version/"

# load inputs ----
load_inputs <- function(){
  SRLI_2024 <<- read_csv("00_raw_data/SRLI_2026.csv")
  ECOrast <<- terra::rast("00_raw_data/ECOrast_p.tif")
  DEMrast <<- terra::rast("00_raw_data/DEM_rast_p.tif")
  ESA_CCI <<- terra::rast("00_raw_data/ESA_CCI_p.tif")
  TDWG_rast <<- terra::rast("00_raw_data/TDWG_rast_p.tif")
  SRLI_dists <<- read_csv("00_raw_data/SRLI_dists.csv")
  #RL_points <<- read_csv("00_raw_data/RL_points.csv")
  #RL_polys <<- sf::st_read("00_raw_data/RL_polys.shp")
  curated_points <<- read_csv("00_raw_data/threatened_occs.csv") %>%
    dplyr::filter(complete.cases(version2))
}

# load parameters ----
load_parameters <- function(){
gbif_key <<- species %>% dplyr::select(gbifKey)
wcvp_id <<- species$WCVP_current_id
print(paste0(wcvp_id))
sp_alt_min <<- species$alt_min
sp_alt_max <<- species$alt_max
tdwg_mask <<- native_mask(ID = wcvp_id, tdwg_level3_sf, SRLI_dists = SRLI_dists)
# tdwg_ids <<-
#   SRLI_dists %>% dplyr::filter(WCVP_current_id == wcvp_id) %>%
#   dplyr::select(area_code_l3) %>%
#   left_join(tdwg_level3, by = c("area_code_l3" = "LEVEL3_COD")) %>%
#   dplyr::select(LEVEL3_NUM)
#tdwg_ids <<- tdwg_ids$LEVEL3_NUM
sp_hab <<- species %>% dplyr::select(ESA_CCI_codes_manual) %>% dplyr::slice_head() %>% as.vector()
sp_hab <<- as.numeric(strsplit(sp_hab$ESA_CCI_codes_manual, ", ")[[1]])
}

world <- ne_countries(scale = "medium", returnclass = "sf")


# get occs ----
get_occs <- function(gbif_key){

  occs <- rgbif::occ_search(
    taxonKey = gbif_key,
    hasGeospatialIssue = FALSE,
    hasCoordinate = TRUE,
    limit = 100000
    # fields = c(
    #   'key',
    #   'scientificName',
    #   'decimalLatitude',
    #   'decimalLongitude',
    #   'coordinateUncertaintyInMeters',
    #   'catalogNumber',
    #   'year'
    # )
  )
  # add coordinate field as sometimes it is not provided
  if (!is.null(occs$data) &
      !("coordinateUncertaintyInMeters" %in% names(occs$data))) {
    occs$data[["coordinateUncertaintyInMeters"]] <- 2000
  }

return(occs$data)

}

# clean occs ----
clean_occs <- function(occs){

  cleaned_occs <- LCr::clean_occs(occs)
  cleaned_occ_sf <- st_as_sf(cleaned_occs, coords =c("decimalLongitude", "decimalLatitude"),
                             crs = getCRS('WGS84'))
  return(cleaned_occ_sf)
}

# native occs ----
native_occ <- function(eco_tdwg, thinned_occs){
  # # global extent, 1 for native range, NA for non-native
  # tdwg <- terra::ifel(test = TDWG_rast %in% tdwg_ids, yes = 1, no = NA)
  #
  # # clip extent to non-na cells
  # tdwg_clip <- clip_rastext(tdwg)

  # get occs that are in native range only
  native_occs <- terra::extract(eco_tdwg, thinned_occs)
  native_occs <- cbind(thinned_occs, native_occs)
  native_occs <- native_occs %>% tidyterra::filter(complete.cases(eco_rast))

  #outputs <- list("native_occs" = native_occs,"tdwg_clip" = tdwg_clip)
  return(native_occs)

}


native_mask <-  function(ID, tdwg_level3, SRLI_dists){

  dist <- SRLI_dists %>%
    dplyr::filter(WCVP_current_id == ID)

  dist_join <- dplyr::left_join(dist, tdwg_level3, by = c("area_code_l3" = "LEVEL3_COD"))

  dist_join = sf::st_sf(geom=dist_join)
}

# ecoregions mask ----
# extract raster values from points
#eco_mask <- function(ECOrast, tdwg_ids, native_occs, tdwg_clip){
eco_mask <- function(eco_tdwg, native_occs){

  eco_ids <- terra::extract(eco_tdwg, native_occs)
  eco_ids <- unique(eco_ids$eco_rast[!is.na(eco_ids$eco_rast)]) # remove NA e.g. sea

  # filter the eco raster
  eco_range <-   terra::ifel(test = eco_tdwg %in% eco_ids,yes = 1,no = NA)

  # crop the eco_range to the TDWG extent
  #eco_crop <- terra::crop(eco_range, tdwg_clip, extend = TRUE)

return(eco_range)

}

# DEM mask
dem_mask <- function(dem_tdwg, native_occs, sp_alt_min, sp_alt_max){

  DEM_range <- dem_tdwg >= sp_alt_min & dem_tdwg <= sp_alt_max

  }


# habitat mask
hab_mask <- function(ESA_tdwg, sp_hab){

  # filter habitat on habitat classes
  HAB_range <-  terra::ifel(test = ESA_tdwg %in% sp_hab, yes = 1, no = NA)

  # rename the rasters to sensible names
  names(HAB_range) <-  paste0("year_", 1992:2020)

  return(HAB_range)

}

# # run_AOH ----
# run_aoh <- function(DEMrast, eco_mask, native_occs, sp_alt_min, sp_alt_max, sp_hab, ESA_CCI){
#
#   # crop DEM
#   #DEM_crop <- terra::crop(DEMrast, eco_mask, mask = TRUE, overwrite = TRUE)
#   #DEM_fltr <- DEM_crop >= sp_alt_min & DEM_crop <= sp_alt_max
#
#   # crop HABITAT
#   HAB_crop <- terra::crop(ESA_CCI, DEM_fltr, overwrite = TRUE)
#   HAB_masked <- terra::mask(HAB_crop, DEM_fltr, maskvalues = TRUE, inverse = TRUE)
#
#   return(HAB_final)
# }



# range_map ----
# range buffer using rapoport
range_map <- function(native_occs, HAB_final){

  # make rapoport buffer, project the points to mollweide
  native_occs_p <- terra::project(native_occs, "+proj=moll")
  native_occs_xy <- as.data.frame(crds(native_occs_p)) %>% dplyr::rename(X = x, Y = y)
  attr(native_occs_xy,'crs') <- "+proj=moll"
  range_buffer <- rapo_mask(native_occs_xy)
  range_buffer <- terra::vect(range_buffer)
  range_clip <- terra::rasterize(range_buffer, HAB_final)

  return(range_clip)
}


# AOH stats ----
# generate statistics e.g. #cell, proportional change and EOO, AOO for raster stack
aoh_stats = function(sp_id, aoh, rl_stats = TRUE) {
  # function to get the cell count per layer
  count_non_na <- function(x)
    sum(!is.na(values(x)))

  # Apply the function to each layer and store results in a vector
  cell_count <- lapply(aoh, count_non_na)

  rast_names <-  paste0("year_", 1992:2020)
  res_df <-
    data.frame(year = rast_names, cell_count = unlist(cell_count))

  # add powo ID
  res_df$powoID <- sp_id

  # add proportion and proportion change
  res_df$proportion <-
    res_df$cell_count / res_df$cell_count[res_df$year == "year_1992"]
  res_df$prop_change <- res_df$proportion - 1

  if (rl_stats == TRUE) {
    # add red list stats
    rl_res <- lapply(aoh, RL_stats)
    rl_res_df <- do.call(dplyr::bind_rows, rl_res)
    res_df <- dplyr::bind_cols(res_df, rl_res_df)
  }

  return(res_df)

}


#[1] "(8/9) make AOH statistics"
#Error in if (EOOArea < 100) { : missing value where TRUE/FALSE needed
#  Calls: aoh_stats -> lapply -> FUN -> <Anonymous>

# RL_stats ----
# Calculate EOO and AOO from AOH
RL_stats = function(aoh_rast) {
  r_eoo_df <- as.data.frame(as.points(aoh_rast), geom = "XY") #aoh_rast
  r_eoo_ll <- r_eoo_df %>%
    dplyr::rename(Y = y) %>%
    dplyr::rename(X = x) %>%
    dplyr::select(X, Y)
  #r_eoo_proj <- rCAT::simProjWiz(r_eoo_ll, returnV = "S")

  # EOO from AOH
  eoo_aoh <- round(rCAT::eoo(r_eoo_ll))
  # catch if EOO is NA
  if (is.na(eoo_aoh)) {
    eoo_aoh_cat = "CR"
  } else {
      eoo_aoh_cat = rCAT::ratingEoo(eoo_aoh)
    }

  # AOO from AOH
  aoo_aoh <- round(rCAT::aoo(r_eoo_ll))
  if (is.na(aoo_aoh)){
    aoo_aoh_cat = "CR"
  } else {
    aoo_aoh_cat = rCAT::ratingAoo(aoo_aoh)
  }

  res_df <- data.frame(
    "eoo_aoh" = eoo_aoh,
    "eoo_aoh_cat" = eoo_aoh_cat,
    "aoo_aoh" = aoo_aoh,
    "aoo_aoh_cat" = aoo_aoh_cat
  )

  return(res_df)

}

# template leaflet map ----
library(leaflet)
library(sf)

leaf_map <-
  function (rast_data = NULL,
            poly_data = NULL,
            point_data = NULL,
            point_lat = NULL,
            point_lng = NULL) {
    map <- leaflet() %>%
      addTiles() %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "esri") %>%
      addScaleBar() %>%
      addLayersControl(
        baseGroups = c("OSM (default)", "esri"),
        overlayGroups = c('poly', 'points', 'raster'),
        options = layersControlOptions(collapsed = FALSE)
      )

    if (!is.null(rast_data)) {
      map <-
        map %>% addRasterImage(rast_data, colors = "green", group = 'raster')
    }

    if (!is.null(poly_data)) {
      map <- map %>% addPolygons(data = poly_data, group = "poly")
    }

    if (!is.null(point_data)) {
      if (inherits(point_data, "sf")) {
        map <-
          map %>% addCircleMarkers(
            data = point_data,
            radius = 3,
            color = "black",
            group = "points"
          )
      } else if (!is.null(point_lat) & !is.null(point_lng)) {
        map <- map %>% addCircleMarkers(
          data = point_data,
          lng = point_lng,
          lat = point_lat,
          radius = 3,
          color = "black",
          group = "points"
        )
      } else {
        stop(
          "Point data must be an sf object or a data frame with specified latitude and longitude columns."
        )
      }
    }

    return(map)
  }

# break raster into chunks ----
# process_chunk <- function(chunk, extent) {
#   cropped_chunk <- crop(chunk, extent)
#   return(cropped_chunk)
# }

# clip raster extent ----
clip_rastext <- function(rast_data) {
  non_na_cells <- which(!is.na(values(rast_data)))
  coords <- terra::xyFromCell(rast_data, non_na_cells)
  xmin <- min(coords[, 1])
  xmax <- max(coords[, 1])
  ymin <- min(coords[, 2])
  ymax <- max(coords[, 2])
  non_na_extent <- terra::ext(xmin, xmax, ymin, ymax)

  cropped_raster <- terra::crop(rast_data, non_na_extent, mask = TRUE)

  return(cropped_raster)

}

# rapoport buffer ----
rapo_mask <- function(thinned_occs) {
  # catch cases when only one point and make 2km buffer
  if (nrow(thinned_occs) <2) {
    occs_spat <- terra::vect(thinned_occs, geom=c("X", "Y"))
    mask_buff <- terra::buffer(occs_spat, 2000)
    mask_buff <- st_as_sf(mask_buff)
  }
  else {
  euMST <- rCAT::eMST(thinned_occs)
  meanMST <- mean(euMST$distance)
  barrierDis <- meanMST * 2
  bufferDis <- meanMST * 2

  # Rapoport
  mask_buff = rCAT::subLocRapoport(
    thinned_occs,
    barrierDis = barrierDis,
    bufferDis = bufferDis,
    returnV = "SF"
  )

  mask_buff <- mask_buff$buffers
  }

  return(mask_buff)
}

# MCP buffer ----
mcp_mask <- function(thinned_occs, buff_dis = NULL) {

  if (nrow(thinned_occs) < 2) {
    # Single point → 2 km buffer
    occs_spat <- terra::vect(thinned_occs, geom=c("X", "Y"))
    mask_buff <- terra::buffer(occs_spat, 2000)
    mask_buff <- st_as_sf(mask_buff)

  } else {
    # Multiple points → MCP
    mcp <- rCAT::eoo(thinned_occs, returnV = "SF")
    mask_buff <- mcp

    # Apply additional user-specified buffer if given
    if (!is.null(buff_dis)) {
      mask_buff <- sf::st_buffer(mask_buff, dist = buff_dis)
    }
  }

  return(mask_buff)
}


prepare_range_input <- function(species, input_data,
                                map_method = c("rapo", "mcp", "mcp_buff"),
                                buff_dis = NULL) {

  # wkt - single or multiple
  if (species$curated_map_source == "wkt") {
    wkt_strings <- na.omit(unique(input_data$errdiam))
    if (length(wkt_strings) == 0) return(NULL)

    range_poly <- vect(wkt_strings, crs = "EPSG:4326")

    # optional: dissolve into one range polygon
    if (length(wkt_strings) > 1) {
      range_poly <- terra::aggregate(range_poly)
    }

    return(range_poly)
  }

  # IUCN points
  if (species$curated_map_source == "iucn_points") {

  # import csv
  iucn_points <- read_csv(paste0("00_raw_data/RL_points/", species$WCVP_current_id, ".csv"))

  # geom
  thinned_occs <- vect(iucn_points,
                       geom = c("dec_long", "dec_lat"),
                       crs = "EPSG:4326")

  thinned_occs <- thinned_occs[thinned_occs$presence == 1, ]

  # project
  native_occs_p <- terra::project(thinned_occs, "+proj=moll")
  native_occs_xy <- as.data.frame(crds(native_occs_p)) %>% dplyr::rename(X = x, Y = y)
  attr(native_occs_xy,'crs') <- "+proj=moll"

  # Select method
  if (map_method == "rapo") {
    range_buffer <- rapo_mask(native_occs_xy)

  } else if (map_method == "mcp") {
    range_buffer <- mcp_mask(native_occs_xy)

  } else if (map_method == "mcp_buff") {
    if (is.null(buff_dis)) {
      stop("For map_method = 'mcp_buff', you must supply buff_dis (in metres).")
    }

    #buff_dis = longestAxis(native_occs_xy) / 10
    #range_buffer <- mcp_mask(native_occs_xy, buff_dis = buff_dis)
    #range_buffer <- mcp_mask(native_occs_xy, buff_dis = longestAxis(native_occs_xy)) / 10
    #print(buff_dis)
    #print(paste0("longest axis = ", longestAxis(native_occs_xy)))

  }

  # Standardise output (always same post-processing)
  occs_poly <- terra::vect(range_buffer)
  crs(occs_poly) <- "+proj=moll"
  range_poly <- terra::project(occs_poly, "EPSG:4326")

  return(range_poly)

  }

  # wkt & occs
  if (species$curated_map_source == "threatened_occs_wkt") {
    # 1. Occurrence points only (exclude WKT rows)
    print("(1/9) Getting occurrences")
    occs <- input_data %>%
      dplyr::filter(is.na(errdiam) | !grepl("POLYGON", errdiam, ignore.case = TRUE))
    #occs <- curated_points %>% dplyr::filter(WCVP_current_id == species$WCVP_current_id)
    occs$decimalLongitude <- occs$longdec
    occs$decimalLatitude <- occs$latdec
    occs$coordinateUncertaintyInMeters <- as.numeric(occs$radius) * 100
    occs <- occs %>% tidyr::unite(scientificName, c(genus, sp1), sep = "_", remove = FALSE)
    occs$wcvp_ipni_id <- species$WCVP_current_id
    occs$GBIF_usageKey <- species$gbifKey

    # bail out if there are no occurrences
    if (nrow(occs) <1) {
      errors <- data.frame(sp_id = wcvp_id, error = "no points")
      write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
      print("no points")
    }  else {

      # step 2 - clean occs
      print("(2/9) Cleaning occurrences")
      #cleaned_occs <- clean_occs(occs)

      native_ranges <- LCr::get_native_range(keys = species)
      checked_occs <- LCr::check_occs(occs, native_ranges, buffer = 2000)
      flagged_occs <- checked_occs$checked_data
      flagged_sum <- checked_occs$summary
      cleaned_occs <- LCr::clean_occs(flagged_occs)
      cleaned_occs <- cleaned_occs$clean_occs

      # bail out if there are no occurrences after cleaning
      if (nrow(cleaned_occs) <1) {
        errors <- data.frame(sp_id = wcvp_id, error = "no curated points after cleaning")
        write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
        print("no points after cleaning")
      }  else {

        # Define required columns and their default type
        required_cols <- c("elevation", "recordNumber", "locality", "habitat")

        # Add missing columns with NA_real_ as default values
        for (col in required_cols) {
          if (!col %in% names(cleaned_occs)) {
            cleaned_occs[[col]] <- NA_real_
          }
        }

        # get the point-based EOO and AOO
        r_eoo_ll <- cleaned_occs %>%
          dplyr::rename(lat = decimalLatitude) %>%
          dplyr::rename(long = decimalLongitude) %>%
          dplyr::select(long, lat)

        r_eoo_proj <- rCAT::simProjWiz(r_eoo_ll, returnV = "S")
        eoo_points <- round(rCAT::eoo(r_eoo_proj))
        aoo_points <- round(rCAT::aoo(r_eoo_proj))
        print(paste0("EOO points: ", eoo_points))
        print(paste0("AOO points: ", aoo_points))

        # elevation check
        if ("elevation" %in% names(cleaned_occs)) {
          elevation <- cleaned_occs %>%
            dplyr::summarise(
              ElevationLower.limit = round(quantile(elevation, 0.05, na.rm = TRUE)),
              ElevationUpper.limit = round(quantile(elevation, 0.95, na.rm = TRUE))
            )
        } else {
          elevation <- data.frame(
            ElevationLower.limit = "",
            ElevationUpper.limit = ""
          )
        }
        print(elevation)

        # bail out if there are no occurrences after cleaning
        if (is.null(cleaned_occs)) {
          errors <- data.frame(sp_id = wcvp_id, error = "no points after cleaning")
          write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
          print("no points after cleaning")
        } else {

          # step 3 - thin occs
          print("(3/9) Thinning occurrences")
          # don't run thin if there are not enough occs (need min 2)
          if (nrow(cleaned_occs) <=2) {
            #thinned_occs <- terra::vect(cleaned_occs)
            thinned_occs <- vect(cleaned_occs,
                                 geom = c("decimalLongitude", "decimalLatitude"),
                                 crs = "EPSG:4326")
          } else {
            # run thin when 2 or more points
            thinned_occs <- vect(cleaned_occs,
                                 geom = c("decimalLongitude", "decimalLatitude"),
                                 crs = "EPSG:4326")

            thinned_occs <- enmSdmX::geoThin(x = thinned_occs, minDist = 2000)
            #thinned_occs <- terra::vect(thinned_occs)
          }

          native_occs_p <- terra::project(thinned_occs, "+proj=moll")
          native_occs_xy <- as.data.frame(crds(native_occs_p)) %>% dplyr::rename(X = x, Y = y)
          attr(native_occs_xy,'crs') <- "+proj=moll"

          # Select method
          if (map_method == "rapo") {
            range_buffer <- rapo_mask(native_occs_xy)

          } else if (map_method == "mcp") {
            range_buffer <- mcp_mask(native_occs_xy)

          } else if (map_method == "mcp_buff") {
            if (is.null(buff_dis)) {
              stop("For map_method = 'mcp_buff', you must supply buff_dis (in metres).")
            }
            range_buffer <- mcp_mask(native_occs_xy, buff_dis = buff_dis)
          }

          # Standardise output (always same post-processing)
          occs_poly <- terra::vect(range_buffer)
          crs(occs_poly) <- "+proj=moll"
          range_poly <- terra::project(occs_poly, "EPSG:4326")

        }
      }
    }

    # 2. WKT polygons only
    wkt_rows <- input_data %>% dplyr::filter(!is.na(errdiam) & grepl("POLYGON", errdiam, ignore.case = TRUE))
    wkt_poly <- NULL
    if (nrow(wkt_rows) > 0) {
      wkt_strings <- unique(wkt_rows$errdiam)
      wkt_poly <- vect(wkt_strings, crs = "EPSG:4326")
      wkt_poly <- terra::aggregate(wkt_poly)  # dissolve overlaps
    }

    # 3. check if polygons exist and have features
    has_occs <- !is.null(occs_poly) && nrow(occs_poly) > 0
    has_wkt  <- !is.null(wkt_poly)  && nrow(wkt_poly)  > 0

    if (has_occs && has_wkt) {
      # make sure both are in the same CRS
      if (crs(occs_poly) != crs(wkt_poly)) {
        occs_poly <- terra::project(occs_poly, crs(wkt_poly))
      }
      # dissolve union into single poly
      range_poly <- terra::union(occs_poly, wkt_poly) |> terra::aggregate()
    } else if (has_occs) {
      range_poly <- terra::aggregate(occs_poly)
    } else if (has_wkt) {
      range_poly <- terra::aggregate(wkt_poly)
    } else {
      range_poly <- NULL
    }
    return(range_poly)

  }

  # IUCN poly
  if (species$curated_map_source == "iucn_polygon") {

    # 1. Locate shapefile
    #   (you may want species-specific paths, e.g. based on wcvp_id or taxon name)
    shp_path <- paste0("00_raw_data/RL_polys/", species$WCVP_current_id, ".shp")

    if (!file.exists(shp_path)) {
      errors <- data.frame(sp_id = wcvp_id, error = "IUCN shapefile not found")
      write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
      print("IUCN shapefile not found")
    } else {

      # 2. Read shapefile into a SpatVector
      iucn_poly <- terra::vect(shp_path)

      # filter on the extant occs
      iucn_poly <- iucn_poly[iucn_poly$PRESENCE == 1, ]

      # 3. Check CRS
      print(crs(iucn_poly))  # should ideally be EPSG:4326 (WGS84)

      # If not WGS84, reproject
      if (is.na(crs(iucn_poly)) || crs(iucn_poly) != "EPSG:4326") {
        iucn_poly <- terra::project(iucn_poly, "EPSG:4326")
      }

      # 4. Dissolve to one polygon (in case multiple features exist)
      range_poly <- terra::aggregate(iucn_poly)

    }
    return(range_poly)
  }

  # threatened_occs
  if (species$curated_map_source == "threatened_occs") {# 1. Occurrence points only (exclude WKT rows)
    print("(1/9) Getting occurrences")
    occs <- input_data %>%
      dplyr::filter(is.na(errdiam) | !grepl("POLYGON", errdiam, ignore.case = TRUE))
    #occs <- curated_points %>% dplyr::filter(WCVP_current_id == species$WCVP_current_id)
    occs$decimalLongitude <- occs$longdec
    occs$decimalLatitude <- occs$latdec
    occs$coordinateUncertaintyInMeters <- as.numeric(occs$radius) * 100
    occs <- occs %>% tidyr::unite(scientificName, c(genus, sp1), sep = "_", remove = FALSE)
    occs$wcvp_ipni_id <- species$WCVP_current_id
    occs$GBIF_usageKey <- species$gbifKey

    # bail out if there are no occurrences
    if (nrow(occs) <1) {
      errors <- data.frame(sp_id = wcvp_id, error = "no points")
      write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
      print("no points")
    }  else {

      # step 2 - clean occs
      print("(2/9) Cleaning occurrences")
      #cleaned_occs <- clean_occs(occs)

      native_ranges <- LCr::get_native_range(keys = species)
      checked_occs <- LCr::check_occs(occs, native_ranges, buffer = 2000)
      flagged_occs <- checked_occs$checked_data
      flagged_sum <- checked_occs$summary
      cleaned_occs <- LCr::clean_occs(flagged_occs)
      cleaned_occs <- cleaned_occs$clean_occs

      # bail out if there are no occurrences after cleaning
      if (nrow(cleaned_occs) <1) {
        errors <- data.frame(sp_id = wcvp_id, error = "no curated points after cleaning")
        write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
        print("no points after cleaning")
      }  else {

        # Define required columns and their default type
        required_cols <- c("elevation", "recordNumber", "locality", "habitat")

        # Add missing columns with NA_real_ as default values
        for (col in required_cols) {
          if (!col %in% names(cleaned_occs)) {
            cleaned_occs[[col]] <- NA_real_
          }
        }

        # get the point-based EOO and AOO
        r_eoo_ll <- cleaned_occs %>%
          dplyr::rename(lat = decimalLatitude) %>%
          dplyr::rename(long = decimalLongitude) %>%
          dplyr::select(long, lat)

        r_eoo_proj <- rCAT::simProjWiz(r_eoo_ll, returnV = "S")
        eoo_points <- round(rCAT::eoo(r_eoo_proj))
        aoo_points <- round(rCAT::aoo(r_eoo_proj))
        print(paste0("EOO points: ", eoo_points))
        print(paste0("AOO points: ", aoo_points))

        # elevation check
        if ("elevation" %in% names(cleaned_occs)) {
          elevation <- cleaned_occs %>%
            dplyr::summarise(
              ElevationLower.limit = round(quantile(elevation, 0.05, na.rm = TRUE)),
              ElevationUpper.limit = round(quantile(elevation, 0.95, na.rm = TRUE))
            )
        } else {
          elevation <- data.frame(
            ElevationLower.limit = "",
            ElevationUpper.limit = ""
          )
        }
        print(elevation)

        # bail out if there are no occurrences after cleaning
        if (is.null(cleaned_occs)) {
          errors <- data.frame(sp_id = wcvp_id, error = "no points after cleaning")
          write_csv(errors, paste0("04_outputs/curated_workflow/errors/errors_", wcvp_id, ".csv"))
          print("no points after cleaning")
        } else {

          # step 3 - thin occs
          print("(3/9) Thinning occurrences")
          # don't run thin if there are not enough occs (need min 2)
          if (nrow(cleaned_occs) <=2) {
            #thinned_occs <- terra::vect(cleaned_occs)
            thinned_occs <- vect(cleaned_occs,
                                 geom = c("decimalLongitude", "decimalLatitude"),
                                 crs = "EPSG:4326")
          } else {
            # run thin when 2 or more points
            thinned_occs <- vect(cleaned_occs,
                                 geom = c("decimalLongitude", "decimalLatitude"),
                                 crs = "EPSG:4326")

            thinned_occs <- enmSdmX::geoThin(x = thinned_occs, minDist = 2000)
            #thinned_occs <- terra::vect(thinned_occs)
          }

          native_occs_p <- terra::project(thinned_occs, "+proj=moll")
          native_occs_xy <- as.data.frame(crds(native_occs_p)) %>% dplyr::rename(X = x, Y = y)
          attr(native_occs_xy,'crs') <- "+proj=moll"

          # # Rapoport buffer method
          # range_buffer <- rapo_mask(native_occs_xy)
          # occs_poly <- terra::vect(range_buffer)
          # crs(occs_poly) <- "+proj=moll"
          # range_poly <- terra::project(occs_poly, "EPSG:4326")
          #
          # # MCP method
          # range_buffer <- mcp_mask(native_occs_xy)
          # occs_poly <- terra::vect(range_buffer)
          # crs(occs_poly) <- "+proj=moll"
          # range_poly <- terra::project(occs_poly, "EPSG:4326")
          #
          # # MCP method with 10km buffer
          # range_buffer <- mcp_mask(native_occs_xy, buff_dis = 10000)
          # occs_poly <- terra::vect(range_buffer)
          # crs(occs_poly) <- "+proj=moll"
          # range_poly <- terra::project(occs_poly, "EPSG:4326")

          # Select method
          if (map_method == "rapo") {
            range_buffer <- rapo_mask(native_occs_xy)

          } else if (map_method == "mcp") {
            range_buffer <- mcp_mask(native_occs_xy)

          } else if (map_method == "mcp_buff") {
            if (is.null(buff_dis)) {
              stop("For map_method = 'mcp_buff', you must supply buff_dis (in metres).")
            }
            range_buffer <- mcp_mask(native_occs_xy, buff_dis = buff_dis)
          }

          # Standardise output (always same post-processing)
          occs_poly <- terra::vect(range_buffer)
          crs(occs_poly) <- "+proj=moll"
          range_poly <- terra::project(occs_poly, "EPSG:4326")

        }
      }
    }

    return(range_poly)

  }

}


# single_validation ----
single_validation <- function(native_occs, HAB_final, dem_tdwg, val_year, wcvp_id) {
  # occurrences minus the validation year
  val_map <- native_occs %>% tidyterra::filter(year != val_year)

  # build rapoport range without validation year occurrences
  native_occs_p <- terra::project(val_map, "+proj=moll")
  native_occs_xy <-  as.data.frame(crds(native_occs_p)) %>% dplyr::rename(X = x, Y = y)
  attr(native_occs_xy, 'crs') <- "+proj=moll"
  val_buffer <- rapo_mask(native_occs_xy)
  val_buffer <- terra::vect(val_buffer)
  # now we need to rasterise the DEM to the shape of the range buffer to get a mask
  #val_clip <- terra::crop(HAB_final$year_1992, val_buffer, mask = TRUE, overwrite = TRUE)
  val_clip <- terra::crop(dem_tdwg, val_buffer, mask = TRUE, overwrite = TRUE)

  # now filter on the correct AOH year to match the validation year
  val_aoh <- HAB_final %>% tidyterra::select(paste0("year_", val_year))
  # filter AOH by range map - this is aoh within the range
  val_aoh <- terra::crop(val_aoh, val_clip, mask = TRUE)

  # see if validation points are intersecting AOH
  # exact, 1km buffer, 2km buffer, 5km buffer
  # filter on the points, project (so that you can build rapo range poly)
  val_yr_occs <- native_occs %>% tidyterra::filter(year == val_year)
  val_yr_occs_p <- terra::project(val_yr_occs, "+proj=moll")
  val_yr_occs_xy <-
    as.data.frame(crds(val_yr_occs_p)) %>% dplyr::rename(X = x, Y = y)
  attr(val_yr_occs_xy, 'crs') <- "+proj=moll"

  # ggplot to illustrate validation
  # ggplot to sense check----
  # print("ggplot for sense check")
  # raster_geog <- terra::project(val_aoh, "EPSG:4326")
  # raster_df <- as.data.frame(raster_geog, xy = TRUE)
  # xmin <- min(raster_df[, 1])
  # xmax <- max(raster_df[, 1])
  # ymin <- min(raster_df[, 2])
  # ymax <- max(raster_df[, 2])
  #
  # val_clip_geog <- terra::project(val_clip, "EPSG:4326")
  # val_clip_df <- as.data.frame(val_clip_geog, xy = TRUE)
  #
  # ggplot() +
  #   geom_sf(data = world) +
  #   geom_tile(data = raster_df, aes(x = x, y = y, fill = "raster_df"), color = "red") +
  #   geom_tile(data = val_clip_df, aes(x = x, y = y, fill = "val_clip_df", alpha = 0.1), color = "blue") +
  #   geom_sf(data = val_yr_occs, color = "black", size = 3, shape = 1) +
  #   coord_sf(
  #     xlim = c(xmin = xmin, xmax = xmax),
  #     ylim = c(ymin = ymin, ymax = ymax)
  #     #xlim = sf::st_bbox(range_clip_geog)[c(1,3)],
  #     #ylim = sf::st_bbox(range_clip_geog)[c(2,4)]
  #     #expand = TRUE
  #   ) +
  #   labs(title = "AOH for validation check") +
  #   theme(legend.position="none")
  # ggsave(paste0("04_outputs/curated_workflow/png_maps/map_", val_year, "_",wcvp_id, ".png"))

  # exact test - do points fall within aoh cells?
  exact_res <- terra::extract(val_aoh, val_yr_occs)

  # 1km buffer test (need to make polygon buffers first)
  buff_1km <- terra::buffer(val_yr_occs, 1000)
  val_res_1km <- terra::extract(val_aoh, buff_1km)
  buff_res_1km <-
    sapply(split(val_res_1km, val_res_1km$ID), function(x) {
      if (any(x[, 2] == 1, na.rm = TRUE)) {
        return(1)
      } else {
        return(NA)
      }
    })
  buff_res_1km_df <-
    data.frame(id = names(buff_res_1km), has_value_1 = buff_res_1km)

  # 2km buffer test
  buff_2km <- terra::buffer(val_yr_occs, 2000)
  val_res_2km <- terra::extract(val_aoh, buff_2km)
  buff_res_2km <-
    sapply(split(val_res_2km, val_res_2km$ID), function(x) {
      if (any(x[, 2] == 1, na.rm = TRUE)) {
        return(1)
      } else {
        return(NA)
      }
    })
  buff_res_2km_df <-
    data.frame(id = names(buff_res_2km), has_value_1 = buff_res_2km)

  # 5km buffer test
  buff_5km <- terra::buffer(val_yr_occs, 5000)
  val_res_5km <- terra::extract(val_aoh, buff_5km)
  buff_res_5km <-
    sapply(split(val_res_5km, val_res_5km$ID), function(x) {
      if (any(x[, 2] == 1, na.rm = TRUE)) {
        return(1)
      } else {
        return(NA)
      }
    })
  buff_res_5km_df <-
    data.frame(id = names(buff_res_5km), has_value_1 = buff_res_5km)

  # summarise results
  val_raw_exact <-
    val_summary(val_year,
                val_clip,
                val_aoh,
                exact_res,
                val_yr_occs,
                extract_method = "exact")
  val_raw_1km <-
    val_summary(val_year,
                val_clip,
                val_aoh,
                buff_res_1km_df,
                val_yr_occs,
                extract_method = "1km_buffer")
  val_raw_2km <-
    val_summary(val_year,
                val_clip,
                val_aoh,
                buff_res_2km_df,
                val_yr_occs,
                extract_method = "2km_buffer")
  val_raw_5km <-
    val_summary(val_year,
                val_clip,
                val_aoh,
                buff_res_5km_df,
                val_yr_occs,
                extract_method = "5km_buffer")

  # final results as df
  raw_val_results <-
    bind_rows(val_raw_exact, val_raw_1km, val_raw_2km, val_raw_5km)

  return(raw_val_results)

}

# validation raw summary ----
val_summary <-
  function(val_year,
           val_clip,
           val_aoh,
           extract_res,
           val_yr_occs,
           extract_method) {
    validation_year <- val_year
    range_pixels <- sum(!is.na(values(val_clip)))
    aoh_pixels <- sum(!is.na(values(val_aoh)))
    aoh_prev <- aoh_pixels / range_pixels
    points_in_aoh <- extract_res %>% dplyr::select(2) %>%
      dplyr::summarise(dplyr::across(everything(), ~ sum(. == 1, na.rm = TRUE))) %>%  dplyr::pull()
    #points_in_aoh <- val_res %>%   select(starts_with("year_")) %>%
    points_out_aoh <- extract_res %>%  dplyr::select(2) %>%
      dplyr::summarise(dplyr::across(everything(), ~ sum(is.na(.)))) %>%  dplyr::pull()
    points_total <- nrow(val_yr_occs)

    # pull it together
    validation_results <- data.frame(
      validation_year = val_year,
      range_pixels = range_pixels,
      aoh_pixels = aoh_pixels,
      aoh_prev = aoh_prev,
      points_in_aoh = points_in_aoh,
      points_out_aoh = points_out_aoh,
      points_total = points_total,
      method = extract_method
    )
    return(validation_results)
  }

# validation test ----
validation_test <- function(validation_scores) {
  validation_scores <- validation_scores %>%
    group_by(method) %>%
    dplyr::summarise(
      meanAOHprev = mean(aoh_prev),
      sumpointsin = sum(points_in_aoh),
      sumpointstot = sum(points_total),
      meanPointsprev = sumpointsin / sumpointstot,
      valid_test =  meanPointsprev > meanAOHprev
    ) %>%
    mutate(sp_id = validation_scores$sp_id[1])

  return(validation_scores)
}



#### extras
# Generate AOH raster stack for 1992 - 2020
aoh_single = function(sp_id,
                      sp_hab,
                      sp_alt_min,
                      sp_alt_max,
                      range_mask,
                      tdwg_level3,
                      DEMrast,
                      ESA_CCI,
                      save_rast = FALSE,
                      save_results = FALSE,
                      rl_stats = FALSE) {
  print(paste0("Processing AOH for ", sp_id))

  # native range mask for final clip
  range_terra <- terra::vect(range_mask$geometry)

  # crop and mask elevation, then crop habitat stack
  print("(1/6) Cropping and masking the DEM")
  DEM_crop <-
    terra::crop(DEMrast, range_terra, mask = TRUE, overwrite = TRUE)

  # filter DEM to upper and lower
  print("(3/6) Filter DEM to min and max elevation")
  DEM_fltr <- DEM_crop >= sp_alt_min & DEM_crop <= sp_alt_max

  #print("(2/6) Cropping the habitat stack")
  HAB_crop <- terra::crop(ESA_CCI, range_terra, overwrite = TRUE)

  # use DEM filter as mask for ESA hab stack
  print("(4/6) Mask habitat stack with elevation")
  HAB_masked <-
    terra::mask(HAB_crop,
                DEM_fltr,
                maskvalues = TRUE,
                inverse = TRUE)

  # final mask
  aoh_masked <- terra::mask(HAB_masked, range_terra)

  # now filter on the habitat classes
  print("(5/6) Filter on the habitat classes")
  aoh <-
    terra::ifel(test = aoh_masked %in% sp_hab,
                yes = 1,
                no = NA)

  names(aoh) <-  paste0("year_", 1992:2020)

  if (save_rast == TRUE) {
    writeRaster(aoh,
                paste0("04_outputs/Eco_ranges/rasters/aoh_", sp_id, ".tif"),
                overwrite = TRUE)

  }

  if (save_results == TRUE) {
    print(paste0("(6/6) Processing stats for: ", sp_id))
    aoh_stats <- aoh_stats(sp_id, aoh, rl_stats = rl_stats)
    #print(paste0("AOH stats= ", aoh_stats))
    write_csv(
      aoh_stats,
      paste0(
        "04_outputs/Eco_ranges/results_csv/results_",
        sp_id,
        ".csv"
      )
    )

  }

}

################
# functions ----
annual_rate <- function(year1, pop1, year2, pop2) {
  yrs <- year2 - year1
  rate <- (pop2 / pop1)^(1 / yrs)
  return(rate)
}

annual_linear_change <- function(year1, pop1, year2, pop2) {
  yrs <- year2 - year1
  change <- (pop2 - pop1) / yrs
  return(change)
}

project_pop_exp <- function(pop_start, rate, years) {
  return(pop_start * rate^years)
}

project_pop_linear <- function(pop_start, change, years) {
  return(pop_start + change * years)
}

assign_category <- function(decline_percent) {
  if (decline_percent <= -80) {
    return("Critically Endangered (CR)")
  } else if (decline_percent <= -50) {
    return("Endangered (EN)")
  } else if (decline_percent <= -30) {
    return("Vulnerable (VU)")
  } else if (decline_percent <= -20) {
    return("Near Threatened (NT)")
  } else if (decline_percent <= -10) {
    return("Possible VU C")
  } else {
    return("Least Concern (LC)")
  }
}

slope_trend <- function(slope, tol = 0.001) {
  if (slope > tol) {
    return("Increasing")
  } else if (slope < -tol) {
    return("Decreasing")
  } else {
    return("Stable")
  }
}

check_model_fit <- function(df, r2_threshold = 0.5) {
  # If all population values are identical, no trend can be estimated
  if (length(unique(df$pop)) == 1) {
    return(list(
      r2_linear = NA_real_,
      r2_exponential = NA_real_,
      slope_linear = 0,
      slope_exponential = 0,
      p_linear = NA_real_,
      p_exponential = NA_real_,
      best_model = "linear",
      reject = TRUE,
      reason = "No change in population"
    ))
  }

  # Fit linear and exponential models
  fit_lin <- lm(pop ~ year, data = df)
  fit_exp <- lm(log(pop) ~ year, data = df)

  r2_lin <- summary(fit_lin)$adj.r.squared
  r2_exp <- summary(fit_exp)$adj.r.squared

  slope_lin <- coef(fit_lin)[2]
  slope_exp <- coef(fit_exp)[2]

  p_lin <- summary(fit_lin)$coefficients["year", "Pr(>|t|)"]
  p_exp <- summary(fit_exp)$coefficients["year", "Pr(>|t|)"]

  # Choose best model based on R², always assign a value
  if (!is.na(r2_exp) && r2_exp >= r2_threshold && r2_exp > r2_lin) {
    best_model <- "exponential"
  } else if (!is.na(r2_lin) && r2_lin >= r2_threshold) {
    best_model <- "linear"
  } else {
    best_model <- "linear" # default if both below threshold
  }

  # Reject if both models below R² threshold
  reject <- (is.na(r2_lin) || r2_lin < r2_threshold) & (is.na(r2_exp) || r2_exp < r2_threshold)

  return(list(
    r2_linear = r2_lin,
    r2_exponential = r2_exp,
    slope_linear = slope_lin,
    slope_exponential = slope_exp,
    p_linear = p_lin,
    p_exponential = p_exp,
    best_model = best_model,
    reject = reject,
    reason = NULL
  ))
}

calc_decline_exp <- function(data, gen_time, assess_year, fit_check_res) {
  period <- 3 * gen_time
  year_past <- assess_year - period

  if (nrow(data) == 2) {
    rate <- annual_rate(data$year[1], data$pop[1], data$year[2], data$pop[2])
    slope <- log(rate)
  } else {
    fit <- lm(log(pop) ~ year, data = data)
    rate <- exp(coef(fit)[2])
    slope <- coef(fit)[2]
  }

  pop_past <- project_pop_exp(data$pop[1], 1 / rate, data$year[1] - year_past)
  pop_present <- project_pop_exp(data$pop[nrow(data)], rate, assess_year - data$year[nrow(data)])

  decline <- (pop_present - pop_past) / pop_past
  decline_percent <- 100 * decline

  return(data.frame(
    method = "exponential",
    generation_time = gen_time,
    assessment_year = assess_year,
    year_past = year_past,
    rate_or_change = rate,
    slope_trend = slope_trend(slope),
    pop_past = pop_past,
    pop_present = pop_present,
    decline_percent = decline_percent,
    category = assign_category(decline_percent),
    r2_linear = fit_check_res$r2_linear,
    r2_exponential = fit_check_res$r2_exponential,
    slope_linear = fit_check_res$slope_linear,
    slope_exponential = fit_check_res$slope_exponential,
    p_linear = fit_check_res$p_linear,
    p_exponential = fit_check_res$p_exponential,
    #aic_best_model = fit_check_res$best_model,
    reject = fit_check_res$reject
  ))
}

calc_decline_linear <- function(data, gen_time, assess_year, fit_check_res) {
  period <- 3 * gen_time
  year_past <- assess_year - period

  if (nrow(data) == 2) {
    change <- annual_linear_change(data$year[1], data$pop[1], data$year[2], data$pop[2])
    slope <- change
  } else {
    fit <- lm(pop ~ year, data = data)
    change <- coef(fit)[2]
    slope <- coef(fit)[2]
  }

  pop_past <- project_pop_linear(data$pop[1], change, year_past - data$year[1])
  pop_present <- project_pop_linear(data$pop[nrow(data)], change, assess_year - data$year[nrow(data)])

  decline <- (pop_present - pop_past) / pop_past
  decline_percent <- 100 * decline

  return(data.frame(
    method = "linear",
    generation_time = gen_time,
    assessment_year = assess_year,
    year_past = year_past,
    rate_or_change = change,
    slope_trend = slope_trend(slope),
    pop_past = pop_past,
    pop_present = pop_present,
    decline_percent = decline_percent,
    category = assign_category(decline_percent),
    r2_linear = fit_check_res$r2_linear,
    r2_exponential = fit_check_res$r2_exponential,
    slope_linear = fit_check_res$slope_linear,
    slope_exponential = fit_check_res$slope_exponential,
    p_linear = fit_check_res$p_linear,
    p_exponential = fit_check_res$p_exponential,
    #aic_best_model = fit_check_res$best_model,
    reject = fit_check_res$reject
  ))
}

plot_decline <- function(data, gen_time, assess_year,species_col = "species",
                         pop_col = "pop") {
  # Run the evaluation once to get results (single species only)
  results <- evaluate_declines(
    data,
    assessment_year = assess_year,
    gen_time = gen_time,
    species_col = species_col,
    pop_col = pop_col
  )

  year_past <- results$year_past
  year_future <- assess_year
  years <- seq(year_past, year_future, by = 1)
  observed_range <- seq(min(data$year), max(data$year), by = 1)

  ## ---- Exponential predictions ----
  if (nrow(data) == 2) {
    rate <- log(data$pop[2]/data$pop[1]) / (data$year[2] - data$year[1])
    pop_exp_all <- project_pop_exp(data$pop[1], exp(rate), years - data$year[1])
    pop_exp_obs <- project_pop_exp(data$pop[1], exp(rate), observed_range - data$year[1])
  } else {
    fit <- lm(log(pop) ~ year, data = data)
    pop_exp_all <- exp(predict(fit, newdata = data.frame(year = years)))
    pop_exp_obs <- exp(predict(fit, newdata = data.frame(year = observed_range)))
  }

  ## ---- Linear predictions ----
  if (nrow(data) == 2) {
    change <- annual_linear_change(data$year[1], data$pop[1], data$year[2], data$pop[2])
    pop_lin_all <- project_pop_linear(data$pop[1], change, years - data$year[1])
    pop_lin_obs <- project_pop_linear(data$pop[1], change, observed_range - data$year[1])
  } else {
    fit <- lm(pop ~ year, data = data)
    pop_lin_all <- predict(fit, newdata = data.frame(year = years))
    pop_lin_obs <- predict(fit, newdata = data.frame(year = observed_range))
  }

  ## ---- Combine prediction data ----
  pred_df <- data.frame(
    year = rep(years, 2),
    pop = c(pop_exp_all, pop_lin_all),
    method = rep(c("Exponential", "Linear"), each = length(years)),
    type = "Predicted"
  )

  obs_fit_df <- data.frame(
    year = rep(observed_range, 2),
    pop = c(pop_exp_obs, pop_lin_obs),
    method = rep(c("Exponential", "Linear"), each = length(observed_range)),
    type = "Observed fit"
  )

  plot_df <- rbind(pred_df, obs_fit_df)

  ## ---- Plot ----
  ggplot() +
    geom_point(data = data, aes(x = year, y = pop), size = 2, color = "black") +
    geom_line(data = plot_df, aes(x = year, y = pop, color = method, linetype = type), size = 1) +
    geom_segment(
      aes(x = year_past, xend = assess_year,
          y = min(data$pop) * 0.95, yend = min(data$pop) * 0.95),
      arrow = arrow(ends = "both", length = unit(0.2, "cm")),
      color = "grey40"
    ) +
    annotate("text",
             x = (year_past + assess_year) / 2,
             y = min(data$pop) * 0.9,
             label = paste0("3-gens (", 3 * gen_time, " yrs)"),
             color = "grey20") +
    labs(
      title = paste0("Observed and Projected Population Trends — ", unique(data[[species_col]])),
      subtitle = paste0("Assessment year: ", assess_year, ", Generation time: ", gen_time),
      x = "Year", y = "Population size", color = "Model", linetype = "Segment"
    ) +
    theme_minimal()
}

evaluate_decline_single <- function(df, assessment_year, gen_time) {

  period <- 3 * gen_time
  year_past <- assessment_year - period

  # constant population
  if(length(unique(df$pop)) == 1){
    return(data.frame(
      species = unique(df$species),
      gen_time = gen_time,
      pop_past = NA_real_,
      year_past = year_past,
      pop_present = NA_real_,
      year_present = assessment_year,
      decline_percent = NA_real_,
      category = NA_character_,
      method = NA_character_,
      slope = 0,
      slope_units = NA_character_,
      slope_trend = "Stable",
      r2 = NA_real_,
      p_value = NA_real_,
      reject = TRUE
    ))
  }

  # two points
  if(nrow(df) == 2){
    slope_lin <- (df$pop[2] - df$pop[1]) / (df$year[2] - df$year[1])
    slope_exp <- log(df$pop[2]/df$pop[1]) / (df$year[2] - df$year[1])

    pop_past <- project_pop_linear(df$pop[1], slope_lin, year_past - df$year[1])
    pop_present <- project_pop_linear(df$pop[2], slope_lin, assessment_year - df$year[2])

    best_model <- "linear"
    slope_val <- slope_lin
    slope_units <- "individuals/year"
    slope_trend_val <- slope_trend(slope_val)
    r2 <- NA
    p_val <- NA
  } else {
    fit_lin <- lm(pop ~ year, data = df)
    slope_lin <- coef(fit_lin)[2]
    r2_lin <- summary(fit_lin)$adj.r.squared
    p_lin <- summary(fit_lin)$coefficients["year","Pr(>|t|)"]
    slope_trend_lin <- slope_trend(slope_lin)

    fit_exp <- lm(log(pop) ~ year, data = df)
    slope_exp <- coef(fit_exp)[2]
    r2_exp <- summary(fit_exp)$adj.r.squared
    p_exp <- summary(fit_exp)$coefficients["year","Pr(>|t|)"]
    slope_trend_exp <- slope_trend(slope_exp)

    if(!is.na(r2_exp) && r2_exp > r2_lin){
      best_model <- "exponential"
      rate <- exp(slope_exp)
      pop_past <- project_pop_exp(df$pop[1], 1/rate, df$year[1] - year_past)
      pop_present <- project_pop_exp(df$pop[nrow(df)], rate, assessment_year - df$year[nrow(df)])
      slope_val <- slope_exp
      slope_units <- "log-growth/year"
      slope_trend_val <- slope_trend_exp
      r2 <- r2_exp
      p_val <- p_exp
    } else {
      best_model <- "linear"
      pop_past <- project_pop_linear(df$pop[1], slope_lin, year_past - df$year[1])
      pop_present <- project_pop_linear(df$pop[nrow(df)], slope_lin, assessment_year - df$year[nrow(df)])
      slope_val <- slope_lin
      slope_units <- "individuals/year"
      slope_trend_val <- slope_trend_lin
      r2 <- r2_lin
      p_val <- p_lin
    }
  }

  decline_percent <- 100 * (pop_present - pop_past) / pop_past
  category <- assign_category(decline_percent)

  data.frame(
    species = unique(df$species),
    gen_time = gen_time,
    pop_past = pop_past,
    year_past = year_past,
    pop_present = pop_present,
    year_present = assessment_year,
    decline_percent = decline_percent,
    category = category,
    method = best_model,
    slope = slope_val,
    slope_units = slope_units,
    slope_trend = slope_trend_val,
    r2 = r2,
    p_value = p_val,
    reject = FALSE
  )
}

evaluate_declines <- function(data,
                              assessment_year = 2015,
                              gen_time = 10,
                              species_col = "species",
                              pop_col = "pop") {

  # Standardize column names once
  df <- data
  colnames(df)[colnames(df) == species_col] <- "species"
  colnames(df)[colnames(df) == pop_col] <- "pop"

  # Ensure vectors
  gen_time_vec <- as.numeric(gen_time)
  assess_year_vec <- as.numeric(assessment_year)

  # Split by species
  species_list <- split(df, df$species)

  results_list <- lapply(species_list, function(subdf) {
    # Loop over assessment years
    res_years <- lapply(assess_year_vec, function(ay) {
      # Loop over generation lengths
      res_gens <- lapply(gen_time_vec, function(gt) {
        evaluate_decline_single(subdf,
                                assessment_year = ay,
                                gen_time = gt)
      })
      do.call(rbind, res_gens)
    })
    do.call(rbind, res_years)
  })

  # Combine all species results
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL
  return(results)
}

