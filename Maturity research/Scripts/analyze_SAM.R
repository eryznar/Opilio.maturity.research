# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
# sdmTMB model
mod <- readRDS("./Maturity research/Data/sdmTMB_spVAR_noBIN_k300.rda")

# SAM
SAM.dat <- read.csv("./Maturity research/Data/SNOW_maleSAM.csv") %>%
  dplyr::rename(SAM = SAM_mean)


ggplot(SAM.dat, aes(YEAR, SAM))+
  geom_line()+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm")

summary(lme(SAM ~ YEAR, data = na.omit(SAM.dat), random = ~ 1 | YEAR, correlation = corAR1()))

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

#Filter predicted specimen data by params (not size yet for full join)
# spec.dat.mat <- spec.dat$specimen %>%
#         filter(YEAR %in% mod$data$YEAR, SHELL_CONDITION == 2, SEX == 1) %>%
#         mutate(SIZE_1MM = floor(SIZE),
#                BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
#                BIN2 = BIN_5MM) %>%
#         separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
#         mutate(LOWER = as.numeric(sub('.', '', LOWER)),
#                UPPER = as.numeric(gsub('.$', '', UPPER)),
#                SIZE_5MM = (UPPER + LOWER)/2,
#                YEAR_SCALED = scale(YEAR)) %>%
#         mutate(SEL = predict(s.gam, newdata= ., type = "response"), # predict size-specific selectivity
#               SAMPLING_FACTOR_SEL = SAMPLING_FACTOR/SEL) %>% # account for size specific selectivity in abundance
#         st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
#         st_transform(., crs = "+proj=utm +zone=2") %>%
#         cbind(st_coordinates(.)) %>%
#         as.data.frame(.) %>%
#         mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
#                LONGITUDE = X/1000,
#                YEAR_F = as.factor(YEAR)) %>%
#         predict(mod, ., type = "response", se = FALSE) %>%
#         rename(PROP_MATURE = est) %>%
#         mutate(SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR_SEL * PROP_MATURE,
#                SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR_SEL-SAMPLING_FACTOR_MATURE)
# 
#  saveRDS(spec.dat.mat, "./Maturity research/Data/sdmTMB_maturespecdat.csv")

spec.dat.mat <- readRDS("./Maturity research/Data/sdmTMB_maturespecdat.csv") # already accounts for selectivity


# Immature abundance
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.dat.mat %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

bioabund.immature.sel <-  crabpack::calc_bioabund(crab_data = mat.dat.sel, species = "SNOW", 
                                        size_min = 40, size_max = 94,  sex = "male", 
                                        shell_condition = c("new_hardshell", "oldshell", "very_oldshell"), years = years) %>%
  group_by(YEAR) %>%
  reframe(IMM_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, IMM_ABUND) 

# All male small crab abundance
bioabund.sm.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                        size_min = 40, size_max = 94,  sex = "male", 
                                        shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(SM_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, SM_ABUND)  %>%
  filter(YEAR >= 1989)

ggplot()+
  geom_line(bioabund.sm.sel, mapping = aes(YEAR, SM_ABUND), color = "black", linewidth = 1)+
  geom_line(bioabund.immature.sel, mapping = aes(YEAR, IMM_ABUND), color = "green", linewidth = 1)+
  theme_bw()

cor(bioabund.sm.sel %>% filter(YEAR %in% years) %>% pull(SM_ABUND), bioabund.immature.sel$IMM_ABUND) # correlated, only using sm abund

# Large male mature abundance
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.dat.mat %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

bioabund.lgmat.sel <-  crabpack::calc_bioabund(crab_data = mat.dat.sel, species = "SNOW", 
                                        size_min = 95, size_max = NULL,  sex = "male", 
                                        shell_condition = c("new_hardshell", "oldshell", "very_oldshell"), years = years) %>%
  group_by(YEAR) %>%
  reframe(MAT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, MAT_ABUND) 


# All male large male abundance
bioabund.lg.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                           size_min = 95, size_max = NULL,  sex = "male", 
                                           shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(LG_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, LG_ABUND) %>%
  filter(YEAR >=1989)

ggplot()+
  geom_line(bioabund.lg.sel, mapping = aes(YEAR, LG_ABUND), color = "black", linewidth = 1)+
  geom_line(bioabund.lgmat.sel, mapping = aes(YEAR, MAT_ABUND), color = "green", linewidth = 1)+
  theme_bw()

cor(bioabund.lg.sel %>% filter(YEAR %in% years) %>% pull(LG_ABUND), bioabund.lgmat.sel$MAT_ABUND) # correlated, only using lg abund


# All male large male abundance
bioabund.all.male.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = NULL, size_max = NULL,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(MALE_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, MALE_ABUND)  %>%
  filter(YEAR >=1989)

# instar 1 abundance (40-60mm) (Sainte Marie?)
instar1 <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 40, size_max = 60,  sex = "male", 
                                            shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(INST1_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, INST1_ABUND)  %>%
  filter(YEAR >= 1989)

abund.dat <- right_join(bioabund.lg.sel, bioabund.sm.sel) %>%
   right_join(., bioabund.all.male.sel) %>%
  right_join(., instar1) %>%
  mutate(PROP_SM = SM_ABUND/MALE_ABUND,
         PROP_LG = LG_ABUND/MALE_ABUND)

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
  arrange(., YEAR) %>%
  dplyr::select(!c(DF_BIOMASS, MALE_ABUND, PROP_LG, PROP_SM, SM_ABUND, X, SPECIES, DISTRICT, SAM_hi, SAM_lo, SAM_sd, VAR_total))

write.csv(model.dat, "./Maturity research/Data/model_dat.csv")

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

# Add in avgs/lags
model.dat2 <- model.dat %>%
              # mutate(INST1_ABUND = scale(INST1_ABUND),
              #        LG_ABUND = scale(LG_ABUND)) %>%
              arrange(., YEAR) %>%
              mutate(INST1_ABUND_AVG2 = zoo::rollmean(INST1_ABUND, k = 2, fill = NA, align = "right"),
                     LG_ABUND_AVG2 = zoo::rollmean(LG_ABUND, k = 2, fill = NA, align = "right"),
                     ICE_AVG2 = zoo::rollmean(ICE, k = 2, fill = NA, align = "right"),
                     # SM_TOCC_AVG2 = zoo::rollmean(SM_TOCC, k = 2, fill = NA, align = "right"),
                     TOCC_AVG2 = zoo::rollmean(TOCC, k = 2, fill = NA, align = "right"),
                     INST1_ABUND_AVG3 = zoo::rollmean(INST1_ABUND, k = 3, fill = NA, align = "right"),
                     LG_ABUND_AVG3 = zoo::rollmean(LG_ABUND, k = 3, fill = NA, align = "right"),
                     ICE_AVG3 = zoo::rollmean(ICE, k = 3, fill = NA, align = "right"),
                     # SM_TOCC_AVG3 = zoo::rollmean(SM_TOCC, k = 3, fill = NA, align = "right"),
                    TOCC_AVG3 = zoo::rollmean(TOCC, k = 3, fill = NA, align = "right"))
                     # #lags
                     # SM_ABUND_LAG2 = lag(SM_ABUND, n = 2),
                     # LG_ABUND_LAG2 = lag(LG_ABUND, n = 2),
                     # ICE_LAG2 = lag(ICE, n = 2),
                     # SM_TOCC_LAG2 = lag(SM_TOCC, n = 2),
                     # LG_TOCC_LAG2 = lag(LG_TOCC, n = 2),
                     # SM_ABUND_LAG3 = lag(SM_ABUND, n =  3),
                     # LG_ABUND_LAG3 = lag(LG_ABUND, n =  3),
                     # ICE_LAG3 = lag(ICE, n =  3),
                     # SM_TOCC_LAG3 = lag(SM_TOCC, n =  3),
                     # LG_TOCC_LAG3 = lag(LG_TOCC, n =  3))



# Fit models
response <- "SAM"

# Generate parameter combos to fit over
lg.pars <-  c(NA, names(model.dat2)[grep("LG_ABUND", names(model.dat2))])
sm.pars <- c(NA, names(model.dat2)[grep("INST1_ABUND", names(model.dat2))])
# sm.tocc.pars <-  c(NA,names(model.dat2)[grep("SM_TOCC", names(model.dat2))])
tocc.pars <-  c(NA,names(model.dat2)[grep("TOCC", names(model.dat2))])
ice.pars <- c(NA, names(model.dat2)[grep("ICE", names(model.dat2))])

combos <- expand_grid(
  #df = df.pars,
  ice = ice.pars,
  lg = lg.pars,
  sm = sm.pars,
  tocc = tocc.pars) %>%
  filter(!(is.na(lg) & is.na(sm) & is.na(ice)))

# Fit models over parameter combinations
fits <- pmap_dfr(
  combos,
  function(lg, sm, ice, tocc) {
    terms <- c(
      if (!is.na(lg)) paste0("s(", lg, ",k=4)") else NULL,
      if (!is.na(sm)) paste0("s(", sm, ",k=4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
      if (!is.na(ice)) paste0("s(", ice, ",k=4)") else NULL
    )
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    fit <- gamm(fml, data = model.dat2, family = gaussian(), 
               method =  "REML", correlation = corAR1())
    tibble(
      sm_term = sm,
      lg_term = lg,
      ice_term = ice,
      tocc_term = tocc,
      k_terms = length(terms),
      AIC = AIC(fit),
      GCV = fit$gcv.ubre,
      edf_total = sum(fit$gam$edf)
    )
  }
)



fits %>%
  arrange(AIC)

mod <- gamm(SAM ~ s(INST1_ABUND, k = 4)+
                 s(LG_ABUND, k = 4)+
                s(TOCC_AVG2, k = 4)+
                 s(ICE_AVG3, k = 4),
                correlation = corAR1(),
                 data = model.dat2,
                 family = gaussian())

diagnose.gamm(mod)
diagnose(mod$gam)


