# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# SAM
SAM.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                  dplyr::select(YEAR, SAM) %>%
                   mutate(sex = "Male"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   filter(YEAR >=1989) %>%
                   mutate(sex = "Female")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, sex = c("Male", "Female")))

m.dat <- SAM.dat %>% filter(sex == "Male")
summary(lme(SAM ~ YEAR, data = na.omit(m.dat), random = ~ 1 | YEAR, correlation = corAR1()))

f.dat <- SAM.dat %>% filter(sex == "Female")
summary(lme(SAM ~ YEAR, data = na.omit(f.dat), random = ~ 1 | YEAR, correlation = corAR1()))

ann_df <- data.frame(
  sex  = c("Male", "Female"),
  x    = c(2000, 2000),   # choose positions you like
  y    = c(75, 40),
  lab  = c("p'=0.02*", "p'=0.88")
)

ggplot(SAM.dat, aes(YEAR, SAM)) +
  geom_line() +
  geom_point() +
  ylab("Size-at-50% maturity (mm)")+
  xlab("Year")+
  geom_smooth(method = "lm") +
  facet_wrap(~ factor(sex, levels = c("Male", "Female")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_SAM_timeseries.png", width = 8, height = 7)


# TARGET_SIZES ----
prop.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_INDPREF) %>%
                   rename(value = PMAT_INDPREF) %>%
                   mutate(type = "Male proportion mature ≥101mm"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_5565) %>%
                   rename(value = PMAT_5565) %>%
                   mutate(type = "Female proportion mature 55-65mm")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, 
                           type = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm")))


m.dat <- prop.dat %>% filter(type == "Male proportion mature ≥101mm")
summary(lme(value ~ YEAR, data = na.omit(m.dat), random = ~ 1 | YEAR, correlation = corAR1()))

f.dat <- prop.dat %>% filter(type =="Female proportion mature 55-65mm")
summary(lme(value ~ YEAR, data = na.omit(f.dat), random = ~ 1 | YEAR, correlation = corAR1()))

ann_df <- data.frame(
  type   = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm"),
  x    = c(2015, 2015),  
  y    = c(0.3, 0.5),
  lab  = c("p'=0.04*", "p'=0.02*")
)

ggplot(prop.dat, aes(YEAR, value)) +
  geom_line() +
  geom_point() +
  ylab("Proportion mature")+
  xlab("Year")+
  geom_smooth(method = "lm") +
  facet_wrap(~ factor(type, levels = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) 

ggsave("./Maturity research/Figures/SNOW_PMAT_timeseries.png", width = 8, height = 7)        

# EBS map
region_layers <- akgfmaps::get_base_layers("sebs")

region_layers$survey.area -> pp
ice <- read.csv("./Maturity research/Output/spatial_ice_means_1980-2025.csv") %>%
  st_as_sf(., coords = c("longitude", "latitude"), crs = crs.latlon) %>%
  st_transform(., st_crs(pp)) %>%
  
  
  ## Load map layers
  map_layers <- readRDS("./Data/map_layers.rda")


## Trim survey grid to survey area
region_layers$survey.grid %>%
  st_transform(crs = st_crs(region_layers$survey.area)) %>%
  st_intersection(region_layers$survey.area) -> survey.grid

## Specify plot boundary, transform to map crs
data.frame(x = c(-178.5, -150), 
           y = c(54.5, 67)) %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(., crs = region_layers$crs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame() -> plot.boundary


## Plot 
ggplot() +
  geom_sf(data = region_layers$bathymetry, color=alpha("grey70")) +
  # geom_sf(data = region_layers$survey.grid, fill=NA, color=alpha("grey70"), linewidth = 1)+
  geom_sf(data = region_layers$survey.area, fill = alpha("cadetblue", alpha=0.3), size = 0) +
  geom_sf(data = region_layers$akland, fill = "grey80", size=0.1) +
  geom_sf(data = region_layers$survey.area, fill = NA) +
  #scale_x_continuous(breaks = c(-180, -175, -170, -165, -160, -155, -150), labels = paste0(c(180, 175, 170, 165, 160, 155, 150), "°W")) + 
  #scale_y_continuous(breaks = c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), labels = paste0(c(52, 54, 56, 58, 60, 62, 64, 66, 68, 70), "°N")) +
  coord_sf(xlim = plot.boundary$X,
           ylim = plot.boundary$Y)+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA, color = "grey70"),
        legend.key.size = unit(0.65,'cm'),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank()) -> study_site

 