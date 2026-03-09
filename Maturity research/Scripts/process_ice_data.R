# PURPOSE: To produce mean Jan-Feb and Mar-Apr sea ice cover for the Bering

# AUTHOR: Emily Ryznar

# 1) LOAD LIBS/PARAMS ----
  
source("./Maturity research/Scripts/load_libs_params.R")

ice.years <- 1980:1988
ice.years <-1989:2000
ice.years <- 2001:2013
ice.years <- 2014:2025

# 2) DOWNLOAD ICE DATA FROM ERA5 (can save each year, but running all years doesn't take long) ----
  # source: https://cds.climate.copernicus.eu/
  
  # specify login credentials for the climate data store
  user_id = "f64421c8-a4c9-4c16-9b25-d2914edc68dc" # this can be found on your user profile
  api_key = "6913841f-d568-40b4-9d22-11da042862f8" # this is the API key, also found on your profile
  
  # set key
  wf_set_key(user = user_id,
             key = api_key) 
  
  # specify request for current year
  request <- list(
    "dataset_short_name" = "reanalysis-era5-single-levels-monthly-means",
    "product_type" = "monthly_averaged_reanalysis",
    "variable" = c("sea_ice_cover"),    
    "year" = ice.years,                     
    "month" = sprintf("%02d", 1:12),
    "day" = sprintf("%02d", 1:31),
    "time" = sprintf("%02d:00", 0:23),
    "area" = c(64, -182, 50, -154),      "format" = "netcdf",                  
    "target" = paste0("ERA5_ice_", min(ice.years), "-", max(ice.years), ".nc") # target file name
  )
  
  # run request (you may need to manually click accept license on website --> follow link in error message if it appears)
  wf_request(
    user     = user_id,
    request  = request,
    transfer = TRUE,
    path     = paste0("./Maturity research/Data/"), # where do you want the data to be saved?
    verbose = TRUE
  )

# 3) PROCESS ICE FILES ----
  # Specify unique ice file names
  files <- list.files("./Maturity research/Data/")
  ice.files <- files[grep("ERA5_ice", files)]
  
  ice.means   <- data.frame()
  ice.spatial <- data.frame()
  
  for(ii in 1:length(ice.files)){
    # Process ice data using tidync()
    tidync(paste0("./Maturity research/Data/", ice.files[ii])) %>%
      hyper_filter(longitude = longitude >= -182 & longitude <= -154,
                   latitude = latitude >= 50 & latitude <= 64) %>%
      activate("siconc") %>%
      hyper_tibble() %>%
      mutate(year = lubridate::year(valid_time),
             month = lubridate::month(valid_time),
             latitude = as.numeric(as.character(latitude)),
             longitude = as.numeric(as.character(longitude))) %>%
      filter(month %in% c(1:4)) -> ice
    
    # non-spatial monthly stats (per year-month)
    ice %>%
      group_by(year, month)  %>%
      reframe(
        n          = dplyr::n(),
        value      = mean(siconc, na.rm = TRUE),
        sd_value   = sd(siconc, na.rm = TRUE),
        se_value   = sd_value / sqrt(n),
        lcl_value  = value - 1.96 * se_value,
        ucl_value  = value + 1.96 * se_value
      ) -> mean.ice
    
    # spatial monthly stats (per year-month-lat-lon)
    ice %>%
      group_by(year, month, latitude, longitude)  %>%
      reframe(
        n          = dplyr::n(),
        value      = mean(siconc, na.rm = TRUE),
        sd_value   = sd(siconc, na.rm = TRUE),
        se_value   = sd_value / sqrt(n),
        lcl_value  = value - 1.96 * se_value,
        ucl_value  = value + 1.96 * se_value
      ) -> spatial.ice
    
    ice.means   <- rbind(ice.means,   mean.ice)
    ice.spatial <- rbind(ice.spatial, spatial.ice)
  }
  
  # compute Jan-Feb and Mar-Apr means (non-spatial) WITH uncertainty

 ice.means <- data.frame()
 ice.spatial <- data.frame()
 for(ii in 1:length(ice.files)){
   # Process ice data using tidync()
   tidync(paste0("./Maturity research/Data/", ice.files[ii])) %>%
     hyper_filter(longitude = longitude >= -182 & longitude <= -154,
                  latitude = latitude >= 50 & latitude <= 64) %>%
     activate("siconc") %>%
     hyper_tibble() %>%
     mutate(year = lubridate::year(valid_time),
            month = lubridate::month(valid_time),
            latitude = as.numeric(as.character(latitude)),
            longitude = as.numeric(as.character(longitude))) %>%
     filter(month %in% c(1:4)) -> ice
   
   ice %>%
     group_by(year, month)  %>%
     reframe(value= mean(siconc)) -> mean.ice
   
   ice %>%
     group_by(year, month, latitude, longitude)  %>%
     reframe(value= mean(siconc)) -> spatial.ice
   
   
   
   ice.means <- rbind(ice.means, mean.ice)
   ice.spatial <- rbind(ice.spatial, spatial.ice)
   
 }

  # Scale, and compute Jan-Feb and Mar-Apr means
  ice.means %>%
    group_by(month) %>%
    mutate(
      name = dplyr::case_when(
        month %in% 1:2 ~ "Jan-Feb ice",
        TRUE          ~ "Mar-Apr ice"
      )
    ) %>%
    ungroup() %>%
    group_by(year, name) %>%
    reframe(
      n          = dplyr::n(),                    # number of monthly values in this seasonal mean
      value      = mean(value, na.rm = TRUE),     # seasonal mean of monthly means
      sd_value   = sd(value, na.rm = TRUE),
      se_value   = sd_value / sqrt(n),
      lcl_value  = value - 1.96 * se_value,
      ucl_value  = value + 1.96 * se_value
    ) -> ice.dat
  
  # compute Jan-Feb and Mar-Apr spatial means WITH uncertainty
  # Scale, and compute Jan-Feb and Mar-Apr means
  ice.spatial %>%
    group_by(month, latitude, longitude) %>%
    mutate(
      name = dplyr::case_when(
        month %in% 1:2 ~ "Jan-Feb ice",
        TRUE          ~ "Mar-Apr ice"
      )
    ) %>%
    ungroup() %>%
    group_by(year, latitude, longitude, name) %>%
    reframe(
      n          = dplyr::n(),
      value      = mean(value, na.rm = TRUE),
      sd_value   = sd(value, na.rm = TRUE),
      se_value   = sd_value / sqrt(n),
      lcl_value  = value - 1.96 * se_value,
      ucl_value  = value + 1.96 * se_value
    ) -> spatial.ice.dat
  
  # Save
  write.csv(ice.dat,
            paste0("./Maturity research/Output/ice_means_1980-", current.year, ".csv"),
            row.names = FALSE)
  
  write.csv(spatial.ice.dat,
            paste0("./Maturity research/Output/spatial_ice_means_1980-", current.year, ".csv"),
            row.names = FALSE)
  
  
  # Get EBS grid for masking
  region_layers <- akgfmaps::get_base_layers("sebs")

  region_layers$survey.area -> pp
  
  ice <- read.csv("./Maturity research/Output/spatial_ice_means_1980-2025.csv") %>%
    st_as_sf(., coords = c("longitude", "latitude"), crs = crs.latlon) %>%
    st_transform(., st_crs(pp)) %>%
    st_intersection(., pp)
  
  # area-mean over EBS WITH uncertainty
  ebs.ice <- ice %>%
    st_drop_geometry() %>%
    na.omit() %>%
    group_by(year, name) %>%
    reframe(
      n          = dplyr::n(),
      value      = mean(value, na.rm = TRUE),
      sd_value   = sd(value, na.rm = TRUE),
      se_value   = sd_value / sqrt(n),
      lcl_value  = value - 1.96 * se_value,
      ucl_value  = value + 1.96 * se_value
    )

  write.csv(ebs.ice,
            paste0("./Maturity research/Output/ebs_ice_means_1980-", current.year, ".csv"),
            row.names = FALSE)
