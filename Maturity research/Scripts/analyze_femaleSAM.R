## ------------------------------------------------------------
# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
# Selectivity ----
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

# Calculate female SAM ----
# Calculate weighted mean SAM for mature female
fem.SAM <- spec.dat.sel$specimen %>% 
              filter(SEX == 2, CLUTCH_SIZE>0) %>%
              group_by(YEAR) %>%
              reframe(SAM = weighted.mean(SIZE, weights = SAMPLING_FACTOR)) %>%
          rbind(., data.frame(YEAR = 2020, SAM = NA))

ggplot(fem.SAM, aes(YEAR, SAM))+
  geom_point()+
  geom_line() +
  theme_bw()


# Weighted prop mature in 55-65
# Calculate weighted mean SAM for mature female
fem.pmat <- spec.dat.sel$specimen %>% 
  filter(SEX == 2, SIZE >=55 & SIZE <=65) %>%
  mutate(MATURE  = case_when(CLUTCH_SIZE >0 ~ 1,
                             TRUE ~ 0)) %>%
  dplyr::select(YEAR, SEX, SIZE, SAMPLING_FACTOR, MATURE) %>%
  group_by(YEAR) %>%
  reframe(TOT_CRAB = sum(SAMPLING_FACTOR),
         MATURE = sum(SAMPLING_FACTOR[MATURE == 1]),
         IMMATURE = TOT_CRAB - MATURE,
         PROP_MATURE = MATURE/TOT_CRAB,
         PROP_IMMATURE = IMMATURE/TOT_CRAB) %>%
  dplyr::select(YEAR, PROP_MATURE) %>%
  rename(PMAT_5565 = PROP_MATURE) 

dat <- right_join(fem.pmat, fem.SAM)

ggplot(dat, aes(YEAR, PMAT_5565))+
  geom_point()+
  geom_line() +
  theme_bw()

summary(lme(PMAT_5565 ~ YEAR, data = na.omit(dat), random = ~ 1 | YEAR, correlation = corAR1()))

# All male large male abundance
bioabund.lg.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 95, size_max = NULL,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(LG_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, LG_ABUND) 

# mature female abundance
bioabund.matfem.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = NULL, size_max = NULL,  sex = "female", 
                                            shell_condition = c("mature_female")) %>%
  group_by(YEAR) %>%
  reframe(MAT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, MAT_ABUND) 


# instar 1 abundance (30-50mm) (Sainte Marie?)
instar1 <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 30, size_max = 50,  sex = "female", 
                                    shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(INST1_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, INST1_ABUND) 

abund.dat <- right_join(bioabund.lg.sel, bioabund.matfem.sel) %>%
  right_join(., instar1) 
unique(is.na(abund.dat))

# industry preferred bioabund for exploitation rate
bioabund.indpref.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                                 size_min = 101, size_max = NULL,  sex = "male") %>%
  mutate(ABUNDANCE = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000) %>% # convert to kt
  dplyr::select(YEAR, ABUNDANCE, BIOMASS)

# Directed fishery data
df.dat <- read.csv("./Maturity research/Data/opilio_directedfishery_catch.csv") %>%
  mutate(directedfish_biomass = Retained_kt+ Discarded_males_kt) %>% #
  dplyr::select(Year, directedfish_biomass) %>%
  rename(YEAR = Year, DF_BIOMASS = directedfish_biomass) %>%
  right_join(., bioabund.indpref.sel) %>% # to calculate exploitation rate
  mutate(DF_BIOMASS = case_when((YEAR %in% c(2020)) ~ NA,
                                (YEAR %in% c(2022:2023)) ~ 0,
                                TRUE ~ DF_BIOMASS),
         EXP_RATE = DF_BIOMASS/BIOMASS) %>%
  dplyr::select(YEAR, DF_BIOMASS, EXP_RATE) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1)))

ggplot(df.dat %>% filter(YEAR >=1990), aes(YEAR, EXP_RATE))+
  geom_point()+
  geom_line()+
  theme_bw()

# Bind with SAM
SAM.abund = right_join(SAM.dat, abund.dat)

unique(is.na(SAM.abund)) # SAM should have NAs in years where no chela were measured
unique(SAM.abund[is.na(SAM.abund$SAM) == TRUE,]$YEAR) #2008, 2012, 2014, 2016

# Load Jan-April ice data
ice <- read.csv(paste0("./Maturity research/Output/ice_means_1989-", current.year, ".csv")) %>%
  #filter(name == "Mar-Apr ice") %>%
  group_by(year) %>%
  reframe(value = mean(value)) %>%
  dplyr::select(year, value) %>%
  rename(YEAR = year, ICE = value)

# Load temperature occupied data
t_occ <- read.csv("./Maturity research/Data/BT_occupied.csv") %>%
  rename(TOCC = temp_occ)

# Bind all dataframes into df for modeling and plot
model.dat <- right_join(SAM.abund, df.dat) %>%
  right_join(., ice) %>%
  right_join(., t_occ %>% dplyr::select(!X)) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1))) %>%
  arrange(YEAR) %>%
  dplyr::select(!c(DF_BIOMASS, MALE_ABUND, PROP_LG, PROP_SM,
                   SM_ABUND, X, SPECIES, DISTRICT,
                   SAM_hi, SAM_lo, SAM_sd, VAR_total))

M <- cor(model.dat %>% dplyr::select(!c(YEAR, SAM)), use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(M,
                   type = "upper",
                   method = "square",
                   order  = "hclust",      # cluster variables
                   addCoef.col = "black") 

mdat.long <- model.dat %>%
  pivot_longer(!YEAR, names_to = "Parameter", values_to = "Value") 

ggplot(mdat.long, aes(YEAR, Value))+
  geom_line()+
  geom_point()+
  facet_wrap(~Parameter, scales = "free_y")+
  theme_bw()

## ------------------------------------------------------------
## 2) Build running means (and keep in one object)
## ------------------------------------------------------------
max_lag <- 6

model.dat2 <- model.dat %>%
  dplyr::select(!EXP_RATE)%>%
  arrange(YEAR) %>%
  mutate(
    # 2‑ and 3‑year running means
    ICE_avg2        = zoo::rollmean(ICE,         k = 2, fill = NA, align = "right"),
    ICE_avg3        = zoo::rollmean(ICE,         k = 3, fill = NA, align = "right"),
    INST1_ABUND_avg2= zoo::rollmean(INST1_ABUND, k = 2, fill = NA, align = "right"),
    INST1_ABUND_avg3= zoo::rollmean(INST1_ABUND, k = 3, fill = NA, align = "right"),
    LG_ABUND_avg2   = zoo::rollmean(LG_ABUND,    k = 2, fill = NA, align = "right"),
    LG_ABUND_avg3   = zoo::rollmean(LG_ABUND,    k = 3, fill = NA, align = "right"),
    TOCC_avg2       = zoo::rollmean(TOCC,        k = 2, fill = NA, align = "right"),
    TOCC_avg3       = zoo::rollmean(TOCC,        k = 3, fill = NA, align = "right")
  )
