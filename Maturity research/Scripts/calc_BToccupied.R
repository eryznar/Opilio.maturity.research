# Calculate "D95" for each size sex group in EBS:
#area of stations that make up 95% of the cumulative cpue

#Author: Erin Fedewa

#Follow ups: 
#move toward spatiotemporal modeling approach to account for changing footprint
#current approach is only using stations sampled every year

# load ----
library(tidyverse)

## Read in setup
#source("./Scripts/get_crab_data.R")
source("./Maturity research/Scripts/load_libs_params.R")

# Selectivity
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  rename(SEL = selectivity, SIZE_5MM = size)

s.gam <- gam(SEL ~ s(SIZE_5MM), data = sel, family = Gamma(link = "log"))

# Specimen data
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")

spec.dat.sel <- spec.dat 
spec.dat.sel$specimen <- spec.dat.sel$specimen %>%
  mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  mutate(SEL = predict(s.gam, newdata= ., type = "response"),
         SAMPLING_FACTOR = SAMPLING_FACTOR/SEL)


# spec.dat$specimen <- spec.dat.mat %>%
#                       dplyr::select(!SAMPLING_FACTOR) %>%
#                       dplyr::rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # I think I want all male temp occupied, not just mature

stations <- read.csv("Y:/KOD_Survey/EBS Shelf/Data_Processing/Data/lookup_tables/station_lookup.csv")

corners <- stations %>% 
  filter(STATION_TYPE == "MTCA_CORNER") %>%
  pull(STATION_ID)


##########################################
# use crabpack to calculate zero-filled, per-station mature CPUE
cpue <- crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                size_min = 40, size_max = NULL,  sex = "male", 
                                shell_condition = c("new_hardshell", "oldshell", "very_oldshell"))  #filtering for relevant range of SAM crab

# cpue.lg <- crabpack::calc_cpue(crab_data = spec.dat, species = "SNOW", 
#                                size_min = 95, size_max = NULL,  sex = "male", 
#                                shell_condition = "new_hardshell")  #filtering for relevant range of SAM crab

#Calculate EBS snow crab temperatures of occupancy (CPUE weighted) 
temp_occ <- cpue%>%
              right_join(., spec.dat$haul) %>%
                dplyr::select(YEAR, LATITUDE, LONGITUDE, STATION_ID, CPUE, GEAR_TEMPERATURE) %>%
                filter(is.na(GEAR_TEMPERATURE) == FALSE) %>%
              group_by(YEAR) %>% 
                summarise(temp_occ = weighted.mean(GEAR_TEMPERATURE, w = CPUE, na.rm = T)) %>%
                filter(is.na(temp_occ) == FALSE) %>%
          filter(YEAR>=1989)

temp_occ %>%
  ggplot(aes(x = YEAR, y = temp_occ))+
  geom_point(size=3)+
  geom_line() +
  theme_bw()




write.csv(temp_occ, "./Maturity research/Data/BT_occupied.csv")
