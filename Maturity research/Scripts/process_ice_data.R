
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
files     <- list.files("./Maturity research/Data/")
ice.files <- files[grep("ERA5_ice", files)]

ice.means   <- data.frame()
ice.spatial <- data.frame()

for (ii in seq_along(ice.files)) {
  tidync(paste0("./Maturity research/Data/", ice.files[ii])) %>%
    hyper_filter(longitude = longitude >= -182 & longitude <= -154,
                 latitude  = latitude  >= 50   & latitude  <= 64) %>%
    activate("siconc") %>%
    hyper_tibble() %>%
    mutate(
      year      = lubridate::year(valid_time),
      month     = lubridate::month(valid_time),
      latitude  = as.numeric(as.character(latitude)),
      longitude = as.numeric(as.character(longitude))
    ) %>%
    filter(month %in% 1:4) -> ice
  
  # monthly means + SD / SE
  ice %>%
    group_by(year, month) %>%
    summarise(
      mean  = mean(siconc, na.rm = TRUE),
      sd    = sd(siconc,   na.rm = TRUE),
      n     = sum(!is.na(siconc)),
      se    = sd / sqrt(n),
      .groups = "drop"
    ) -> mean.ice
  
  # cell-level monthly means
  ice %>%
    group_by(year, month, latitude, longitude) %>%
    summarise(
      mean = mean(siconc, na.rm = TRUE),
      .groups = "drop"
    ) -> spatial.ice
  
  ice.means   <- bind_rows(ice.means,   mean.ice)
  ice.spatial <- bind_rows(ice.spatial, spatial.ice)
}

## Jan–Feb / Mar–Apr means with SE (per grid cell and per year)
ice.dat <- ice.means %>%
  mutate(name = if_else(month %in% 1:2, "Jan-Feb ice", "Mar-Apr ice")) %>%
  group_by(year, name) %>%
  summarise(
    value = mean(mean),
    se    = sqrt(sum(se^2) / n()),
    .groups = "drop"
  )

spatial.ice.dat <- ice.spatial %>%
  mutate(name = if_else(month %in% 1:2, "Jan-Feb ice", "Mar-Apr ice")) %>%
  group_by(year, latitude, longitude, name) %>%
  summarise(
    value = mean(mean),
    .groups = "drop"
  )

write.csv(
  ice.dat,
  paste0("./Maturity research/Output/ice_means_1980-", current.year, ".csv"),
  row.names = FALSE
)
write.csv(
  spatial.ice.dat,
  paste0("./Maturity research/Output/spatial_ice_means_1980-", current.year, ".csv"),
  row.names = FALSE
)

## EBS‑wide annual means with SE directly from spatial.ice.dat
region_layers$survey.area -> pp

ice_sf <- spatial.ice.dat %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = crs.latlon) %>%
  st_transform(st_crs(pp)) %>%
  st_intersection(pp)

ebs.ice <- ice_sf %>%
  as.data.frame() %>% 
  group_by(year) %>%          # keep Jan–Feb vs Mar–Apr if you want
  summarise(
    value = mean(value, na.rm = TRUE),
    sd    = sd(value,   na.rm = TRUE),
    n     = sum(!is.na(value)),
    se    = sd / sqrt(n),
    .groups = "drop"
  )

write.csv(
  ebs.ice,
  paste0("./Maturity research/Output/ebs_ice_means_1980-", current.year, ".csv"),
  row.names = FALSE
)