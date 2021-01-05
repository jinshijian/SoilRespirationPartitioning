# load packages
library(tibble)
library(lubridate)
library(kableExtra)
library(raster)
library(ncdf4)
library(tidyr)
library(dplyr)
library(magrittr)
library(data.table)

DATA_DIR <- 'Data'

#*****************************************************************************************************************
# plan 
#*****************************************************************************************************************
plan <- drake::drake_plan(
  srdb = read.csv(here::here("Data", "srdb-data.csv")) %>%
    mutate(RC_annual2 = Ra_annual / Rs_annual,
           RC_annual3 = (Rs_annual - Rh_annual) / Rs_annual) %>% 
    mutate(RC_annual = coalesce(RC_annual, RC_annual2, RC_annual3)) %>% 
    dplyr::select(Site_ID, Latitude, Longitude, Leaf_habit, MAT, MAP, Ecosystem_type, # ecosystem_type reported in the paper
                  Rs_annual, Rh_annual, RC_annual, Manipulation,
                  Partition_method) ,
    # filter(RC_annual >= 0, RC_annual <= 1, Manipulation == "None"
    #        !is.na(Site_ID), !is.na(Latitude), !is.na(Longitude),
    #        !is.na(Leaf_habit), !is.na(Rs_annual), !is.na(Rh_annual)
    #        ),
  
  # Extract Climate region
  IGBP_Koppen_MODIS = read.csv(here::here("Data", "IGBP_Koppen_MODIS.csv")) %>% 
    # Regroup climate data into fewer categories for easier analysis
    mutate(MiddleClimate = case_when(
      ClimateTypes %in% c("Af", "Am", "As", "Aw") ~ "A",
      ClimateTypes %in% c("BSh", "BSk", "BWh", "BWk") ~ "B",
      ClimateTypes %in% c("Cfa", "Cfb", "Cfc") ~ "Cf",
      ClimateTypes %in% c("Csa", "Csb", "Csc") ~ "Cs",
      ClimateTypes %in% c("Cwa", "Cwb", "Cwc") ~ "Cw",
      ClimateTypes %in% c("Dfa", "Dfb", "Dfc", "Dfd") ~ "Df",
      ClimateTypes %in% c("Dsa", "Dsb", "Dsc", "Dwa", "Dwb", "Dwc", "Dwd") ~ "Dsw",
      ClimateTypes %in% c("EF", "ET") ~ "E",
      TRUE ~ "Other")), 
  
  precip = getData("worldclim", path = here::here(), var = "prec", res = 10, download = !file.exists("wc10/prec1.hdr")),
  precip_agg = aggregate(precip, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  tmean = getData("worldclim", path = here::here(), var = "tmean", res = 10, download = !file.exists("wc10/wc10/tmean1.hdr")),
  tmean_agg = aggregate(tmean, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  
  # Calculate annual sum for precipitation...
  precip_global = raster::as.data.frame(precip, xy = TRUE) %>% drop_na(),
  map_global = precip_global %>% dplyr::select(-x, -y) %>% apply(1, sum),
  
  # Calculate annual sum for temperature...
  tmean_global = raster::as.data.frame(tmean, xy = TRUE) %>% drop_na(), 
  mat_global = tmean_global %>% dplyr::select(-x, -y) %>% apply(1, mean), 
  
  mat = tibble(x = tmean_global$x, y = tmean_global$y, mat = as.vector(mat_global)), 
  map = tibble(x = precip_global$x, y = precip_global$y, map = as.vector(map_global)),  
  
  map_mat_global = left_join(map, mat, by = c("x", "y")), 
  
  # MAP data that matches the srdb coordinates
  precip_coords = raster::extract(precip, srdb %>% dplyr::select(Longitude, Latitude)),
  precip_coords_global = raster::extract(precip_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)), # aggregate to 0.5 resolution
  MAP_WC = apply(precip_coords, 1, sum),
  MAP_WC_global = apply(precip_coords_global, 1, sum), 
  
  # The same for MAT
  tmean_vals = raster::extract(tmean, srdb %>% dplyr::select(Longitude, Latitude)),
  tm_coords_global = raster::extract(tmean_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)), # aggregate to 0.5 resolution
  MAT_WC = apply(tmean_vals, 1, mean),
  MAT_WC_global = apply(tm_coords_global, 1, mean), 
  
  ## Pulling Mycorrhizae Data and Investigating it
  am = raster(here::here("Data", "MycDistrAM_current.TIF")),
  am_agg = aggregate(am, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  em = raster(here::here("Data", "MycDistrEM_current.TIF")),
  em_agg = aggregate(em, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  er = raster(here::here("Data", "MycDistrER_current.TIF")),
  er_agg = aggregate(er, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  nm = raster(here::here("Data", "MycDistrNM_current.TIF")),
  nm_agg = aggregate(nm, fact = c(3, 3), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  
  # Extract the myco data for each coordinate pair in the SRDB database and add it to the rest of the data
  AM_percent = raster::extract(am, srdb %>% dplyr::select(Longitude, Latitude)),
  EM_percent = raster::extract(em, srdb %>% dplyr::select(Longitude, Latitude)),
  ER_percent = raster::extract(er, srdb %>% dplyr::select(Longitude, Latitude)),
  NM_percent = raster::extract(nm, srdb %>% dplyr::select(Longitude, Latitude)),
  
  AM_percent_global = raster::extract(am_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  EM_percent_global = raster::extract(em_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  ER_percent_global = raster::extract(er_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  NM_percent_global = raster::extract(nm_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # Extract biomass
  BM_aboveground = raster(here::here("Data", "aboveground_biomass_carbon_2010.tif")) %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  BM_belowground = raster(here::here("Data", "belowground_biomass_carbon_2010.tif")) %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  
  BMa = raster(here::here("Data", "aboveground_biomass_carbon_2010.tif")),
  BMa_agg = aggregate(BMa, fact = c(180, 180), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  BMb = raster(here::here("Data", "belowground_biomass_carbon_2010.tif")),
  BMb_agg = aggregate(BMb, fact = c(180, 180), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  BM_aboveground_global = raster::extract(BMa_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  BM_belowground_global = raster::extract(BMb_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # Add these biomass values to the larger dataset
  # Remove values where::here there::here are no biomass data
  
  ## Extracts N deposition data for each srdb point, already in 0.5 resolution
  #  https://www.isimip.org/gettingstarted/availability-input-data-isimip2b/
  N_dep_half_deg = raster(here::here("Data", "global_mean_Nitrogen_depostion_1980_2017_half_degree.tif")) %>% 
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  
  N_dep_half_deg_global = raster(here::here("Data", "global_mean_Nitrogen_depostion_1980_2017_half_degree.tif")) %>% 
    raster::extract(IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # Pull N-dep values for every coordinate pair in the SRDB
  N_dep_1993 = raster(here::here("Data", "sdat_830_2_20200721_153826639.asc")) %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  
  N_dep_1993_global = raster(here::here("Data", "sdat_830_2_20200721_153826639.asc")) %>%
    raster::extract(IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  ## Extracts GPP from FLUXCOM (already 0.5 resolution)
  gpp_fluxcom = raster(here::here("Data", "fluxcom_gpp2.tif")) %>% 
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude )),
  
  gpp_fluxcom_global = raster(here::here("Data", "fluxcom_gpp2.tif")) %>% 
    raster::extract(IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # Add this to the data pool and remove entries without N-dep data
  srdb_NDep = cbind(srdb, MAP_WC, MAT_WC, AM_percent, EM_percent, ER_percent, NM_percent, BM_aboveground, BM_belowground,
                    N_dep_1993, N_dep_half_deg, gpp_fluxcom),
  
  # Change Latitude and Longitude to the same 0.5*0.5 resolution as in the dataset
  srdb_IGBP = srdb_NDep %>% 
    mutate(Latitude2 = round(Latitude * 2) / 2 + 0.25, 
           Longitude2 = round(Longitude * 2) / 2 + 0.25) %>% 
    # Add data to the large dataset
    left_join(IGBP_Koppen_MODIS, by = c("Latitude2" = "Latitude",
                                        "Longitude2" = "Longitude")) %>% 
    dplyr::select(-Latitude2, -Longitude2), 
  
  # extract BD, clay percentage, and SOC from SoilGrids
  bd_soilgrids = raster("Data/BLDFIE_M_sl2_1km_ll.tif") %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  clay_soilgrids = raster("Data/CLYPPT_M_sl2_1km_ll.tif") %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  soc_soilgrids = raster("Data/OCSTHA_M_sd2_1km_ll.tif") %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  
  soilbd = raster("Data/BLDFIE_M_sl2_1km_ll.tif"),
  soilbd_agg = aggregate(soilbd, fact = c(60, 60), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  bd_soilgrids_global = raster::extract(soilbd_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  soilclay = raster("Data/CLYPPT_M_sl2_1km_ll.tif"),
  soilclay_agg = aggregate(soilclay, fact = c(60, 60), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  clay_soilgrids_global = raster::extract(soilclay_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  soilsoc = raster("Data/OCSTHA_M_sd2_1km_ll.tif"),
  soilsoc_agg = aggregate(soilsoc, fact = c(60, 60), fun = mean, na.rm = TRUE), # aggregate to 0.5 resolution
  soc_soilgrids_global = raster::extract(soilsoc_agg, IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # get EVI 
  EVI_mean = raster("Data/EVI_mean.tif") %>%
    raster::extract(srdb %>% dplyr::select(Longitude, Latitude)),
  
  EVI_mean_global = raster("Data/EVI_mean.tif") %>%
    aggregate(fact = c(5, 5), fun = mean, na.rm = TRUE) %>% # aggregate to 0.5 resolution
    raster::extract(IGBP_Koppen_MODIS %>% dplyr::select(Longitude, Latitude)),
  
  # srdb with all factors extracted
  srdb_all = cbind(srdb_IGBP, bd_soilgrids, clay_soilgrids, soc_soilgrids, EVI_mean) %>% 
    mutate(bd_soilgrids = bd_soilgrids / 1000, # unit from kg/m3 to g/cm3 
           MAT_WC = MAT_WC / 10), # Temp data is stored in degC * 10, so we need to divide 10 to get back to degC
  
  # prepare global 0.5*0.5 resolution data for global RC_prediction map
  globalData = cbind(IGBP_Koppen_MODIS, MAP_WC_global, MAT_WC_global, AM_percent_global, EM_percent_global, 
                     ER_percent_global, NM_percent_global, BM_aboveground_global, BM_belowground_global,
                     N_dep_1993_global, N_dep_half_deg_global, gpp_fluxcom_global,
                     bd_soilgrids_global, clay_soilgrids_global, soc_soilgrids_global, EVI_mean_global) %>% 
    dplyr::rename(MAT = MAT_WC_global, MAP = MAP_WC_global,
                  AM_percent = AM_percent_global, EM_percent = EM_percent_global,
                  ER_percent = ER_percent_global, NM_percent = NM_percent_global, 
                  BM_aboveground = BM_aboveground_global, BM_belowground = BM_belowground_global,
                  N_dep_1993 = N_dep_1993_global, N_dep_half_deg = N_dep_half_deg_global, gpp_fluxcom = gpp_fluxcom_global, 
                  bd_soilgrids = bd_soilgrids_global, clay_soilgrids = clay_soilgrids_global,
                  soc_soilgrids = soc_soilgrids_global,
                  EVI_mean = EVI_mean_global) %>% 
    dplyr::rename(Rs_warner = warner_rs) %>% 
    mutate(MiddleClimate = as.factor(MiddleClimate),
           bd_soilgrids = bd_soilgrids / 1000, # unit from kg/m3 to g/cm3 
           MAT = MAT / 10), # Temp data is stored in degC * 10, so we need to divide 10 to get back to degC
  
  # Global map
  counties = map_data("world", region = ".", exact = FALSE)
)

drake::make(plan)
