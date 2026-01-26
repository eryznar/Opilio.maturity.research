# PURPOSE: To produce mean Jan-Feb and Mar-Apr sea ice cover for the Bering

# AUTHOR: Emily Ryznar

# 1) LOAD LIBS/PARAMS ----
  
source("./Maturity data processing/Scripts/1) load_libs_params.R")

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
    "area" = c(64, -180, 55, -160),           # Bering 
    "format" = "netcdf",                  
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
 
 ice.means <- data.frame()
 for(ii in 1:length(ice.files)){
   # Process ice data using tidync()
   tidync(paste0("./Maturity research/Data/", ice.files[ii])) %>%
     hyper_filter(longitude = longitude >= -180 & longitude <= -165,
                  latitude = latitude >= 55 & latitude <= 63) %>%
     activate("siconc") %>%
     hyper_tibble() %>%
     mutate(year = lubridate::year(valid_time),
            month = lubridate::month(valid_time),
            latitude = as.numeric(as.character(latitude)),
            longitude = as.numeric(as.character(longitude))) %>%
     filter(month %in% c(1:4)) %>% #months = Jan-Apr
     group_by(year, month)  %>%
     reframe(value= mean(siconc)) -> ice
   
   
   ice.means <- rbind(ice.means, ice)
   
 }

  # Scale, and compute Jan-Feb and Mar-Apr means
  ice.means %>%
    group_by(month) %>%
    mutate(value = scale(value),
           name = case_when((month %in% 1:2) ~ "Jan-Feb ice",
                            TRUE ~ "Mar-Apr ice")) %>%
    ungroup() %>%
    group_by(year, name) %>%
    reframe(value = mean(value)) -> ice.dat
    

  # Save
  write.csv(ice.dat, paste0("./Maturity research/Output/ice_means_1989-", current.year, ".csv"), row.names = FALSE)
  