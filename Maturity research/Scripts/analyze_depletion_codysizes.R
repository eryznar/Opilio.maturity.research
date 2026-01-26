# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:
# - Prop mature by small and large bins using chela data? or specimen? or biomass/abundance?


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/1) load_libs_params.R")

source("./Maturity research/Scripts/mature_cpue_bioabund_function.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)


# Specimen data
spec.dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")

# industry preferred bioabund for exploitation rate
bioabund.indpref <-  crabpack::calc_bioabund(crab_data = spec.dat, species = "SNOW", 
                                             size_min = 101, size_max = NULL,  sex = "male") %>%
  mutate(ABUNDANCE = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000) %>% # convert to kt
  dplyr::select(YEAR, ABUNDANCE, BIOMASS)

# Directed fishery data
df.dat <- read.csv("./Maturity research/Data/opilio_directedfishery_catch.csv") %>%
  mutate(directedfish_biomass = Retained_kt+ Discarded_males_kt) %>% #
  dplyr::select(Year, directedfish_biomass) %>%
  rename(YEAR = Year, DF_BIOMASS = directedfish_biomass) %>%
  right_join(., bioabund.indpref) %>% # to calculate exploitation rate
  mutate(DF_BIOMASS = case_when((YEAR %in% c(2020, 2022:2023)) ~ NA,
                                TRUE ~ DF_BIOMASS),
         EXP_RATE = DF_BIOMASS/BIOMASS) %>%
  dplyr::select(YEAR, DF_BIOMASS, EXP_RATE)

ggplot(df.dat, aes(YEAR, EXP_RATE))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_smooth(method = "lm")

# sdmTMB model
mod <- readRDS(paste0(remote_dir, "./SNOW/sdmTMB/s(SIZE)_iid_200_sdmTMB.rda")) 

# Assess size bins
ogive.dat <- rbind(legacy.ogives, sdmTMB.ogives, gam.ogives)

ggplot(ogive.dat %>% filter(Estimator == "sdmTMB"), aes(SIZE, PROP_MATURE, group = YEAR, color = YEAR))+
  geom_line(linewidth = 1)+
  # scale_color_manual(values = c(
  #   "sdmTMB"       = "cadetblue"), name = "")+
  # scale_fill_manual(values = c(
  #   "sdmTMB"       = "cadetblue"), name = "")+
  #facet_wrap(~YEAR)+
  geom_vline(xintercept = 55, linetype = "dashed")+
  geom_vline(xintercept = 65, linetype = "dashed")+
  geom_vline(xintercept = 75, linetype = "dashed")+
  geom_vline(xintercept = 95, linetype = "dashed")+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 175, by = 10))+
  xlim(c(25, 140))+
  ylab("Proportion mature")+
  xlab("Carapace width (mm)")+
  theme(
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12))

ggsave("./Maturity research/Figures/ogive_comparison_sdmTMB_binexplore.png", width  = 10, height = 6, units = "in")


# Filter predicted specimen data by params (not size yet for full join)
preds <- spec.dat$specimen %>%
  filter(YEAR %in% years, SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2,
         BIN = case_when((SIZE_1MM<=95) ~ "Small", # apply small and large bins for modeling below
                         (SIZE_1MM>95) ~ "Large")) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR)) %>%
  predict(mod, ., type = "response", se = FALSE) %>%
  rename(PROP_MATURE = est) %>%
  mutate(SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR * PROP_MATURE)


mod.dat <- preds %>%
  filter(is.na(BIN) == FALSE) %>% # filtering out crab not in those bins) 
  rename(MAT_COUNT = SAMPLING_FACTOR_MATURE,
         COUNT = SAMPLING_FACTOR) %>%
  group_by(YEAR, BIN) %>%
  reframe(COUNT = sum(COUNT),
          MAT_COUNT = sum(MAT_COUNT),
          PROP_MATURE = MAT_COUNT/COUNT) 

mod.dat2 <- mod.dat %>% 
  right_join(., data.frame(YEAR = rep(seq(min(.$YEAR), max(.$YEAR)), 2), BIN = c("Large", "Small")))

ggplot(mod.dat2, aes(YEAR, PROP_MATURE, color = BIN, fill= BIN))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  scale_fill_manual(values = c("navy", "sienna"), labels = c("Large (95-115mm)", "Small (55-75mm)"))+
  scale_color_manual(values = c("navy", "sienna"), labels = c("Large (95-115mm)", "Small (55-75mm)"))+
  theme_bw()+
  ggtitle("Proportion mature by bin")+
  geom_smooth(method = "lm", alpha = 0.15)

# lg.dat$YEAR <- as.numeric(as.character(lg.dat$YEAR))
m1 <- lme(PROP_MATURE ~ YEAR, random = ~ 1 | YEAR, data = mod.dat %>% filter(BIN == "Large"), correlation = corAR1())
m2 <-  lme(PROP_MATURE ~ YEAR, random = ~ 1 | YEAR, data = mod.dat %>% filter(BIN == "Small"), correlation = corAR1())

# Get mature biomass and abundance for large and small bins
ba.dat <- spec.dat
ba.dat$specimen <- preds %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # renaming mature SF to SF so crabpack recognizes

bioabund.sm <-  crabpack::calc_bioabund(crab_data = ba.dat, species = "SNOW", years = years, 
                                        size_min = NULL, size_max = 95,  sex = "male", 
                                        shell_condition = "new_hardshell") %>%
  mutate(MAT_ABUND = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000) %>% # convert to kt
  rename(MAT_ABUND_SM = MAT_ABUND, BIOMASS_SMALL  = BIOMASS) %>%
  dplyr::select(YEAR, MAT_ABUND_SM, BIOMASS_SMALL)

bioabund.lg <-  crabpack::calc_bioabund(crab_data = ba.dat, species = "SNOW", years = years, 
                                        size_min = 95, size_max = NULL,  sex = "male", 
                                        shell_condition = "new_hardshell") %>%
  mutate(MAT_ABUND = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000) %>% # convert to kt
  rename(MAT_ABUND_LG = MAT_ABUND, BIOMASS_LARGE  = BIOMASS) %>%
  dplyr::select(YEAR, MAT_ABUND_LG, BIOMASS_LARGE)

abund.bin.dat <- right_join(bioabund.sm, bioabund.lg) %>%
  right_join(mod.dat, .) 

# Load March-April ice data
ice <- read.csv(paste0("./Maturity research/Output/ice_means_1989-", current.year, ".csv")) %>%
  #filter(name == "Mar-Apr ice") %>%
  group_by(year) %>%
  reframe(value = mean(value)) %>%
  dplyr::select(year, value) %>%
  rename(YEAR = year, ICE = value)

# SAM
SAM <- read.csv("./Maturity research/Output/SAM_comparison.csv") %>%
  filter(SAM_type == "ogive_SAM") %>%
  dplyr::select(!c(SAM_type, X))


# Bind all dataframes into df for modeling
model.dat <- right_join(df.dat, abund.bin.dat) %>%
  right_join(., ice) %>%
  right_join(., SAM) %>%
  group_by(BIN) %>%
  mutate(MAT_ABUND_SM_LAG1 = lag(MAT_ABUND_SM, n=1),
         MAT_ABUND_SM_LAG2 = lag(MAT_ABUND_SM, n=2),
         MAT_ABUND_LG_LAG1 = lag(MAT_ABUND_LG, n=1),
         MAT_ABUND_LG_LAG2 = lag(MAT_ABUND_LG, n=2),
         EXP_RATE_LAG1 = lag(EXP_RATE, n=1),
         EXP_RATE_LAG2 = lag(EXP_RATE, n=2),
         ICE_LAG3 = lag(ICE, n=3),
         #MAT_ABUND_SM_AVG1 = MAT_ABUND_SM,
         MAT_ABUND_SM_AVG2 = zoo::rollmean(MAT_ABUND_SM, k = 2, fill = NA, align = "center"),
         MAT_ABUND_SM_AVG3 = zoo::rollmean(MAT_ABUND_SM, k = 3, fill = NA, align = "center"),
         #MAT_ABUND_LG_AVG1 = MAT_ABUND_LG,
         MAT_ABUND_LG_AVG2 = zoo::rollmean(MAT_ABUND_LG, k = 2, fill = NA, align = "center"),
         MAT_ABUND_LG_AVG3 = zoo::rollmean(MAT_ABUND_LG, k = 3, fill = NA, align = "center"),
         #EXP_RATE_AVG1 = EXP_RATE,
         EXP_RATE_AVG2 = zoo::rollmean(EXP_RATE, k = 2, fill = NA, align = "center"),
         EXP_RATE_AVG3 = zoo::rollmean(EXP_RATE, k = 3, fill = NA, align = "center"),
         #ICE_AVG1 = ICE,
         ICE_AVG2 = zoo::rollmean(ICE, k = 2, fill = NA, align = "center"),
         ICE_AVG3 = zoo::rollmean(ICE, k = 3, fill = NA, align = "center")) #%>%
# rename(MAT_ABUND_SM_LAG0 = MAT_ABUND_SM,
#        MAT_ABUND_LG_LAG0 = MAT_ABUND_LG,
#        EXP_RATE_LAG0 = EXP_RATE, # check dates, usually January before the summer survey
#        ICE_LAG0 = ICE)


#write.csv(model.dat, "./Maturity research/Data/snow_male_GAM_modeldat.csv")

# PROP_MATURE in SMALL/LARGE BIN ----
# Filter dat
sm.dat <- model.dat %>% filter(BIN == "Small")
lg.dat <- model.dat %>% filter(BIN == "Large")

response <- "PROP_MATURE"

# Lags and smooths
abund.lg.pars <- c(c(NA, names(model.dat)[grep("MAT_ABUND_LG_LAG", names(model.dat))]), # lags
                   c(names(model.dat)[grep("MAT_ABUND_LG_AVG", names(model.dat))])) # smooths
abund.sm.pars <- c(c(NA, names(model.dat)[grep("MAT_ABUND_SM_LAG", names(model.dat))]), # lags
                   c(names(model.dat)[grep("MAT_ABUND_SM_AVG", names(model.dat))])) # smooths
df.pars <- c(c(NA,names(model.dat)[grep("EXP_RATE_LAG", names(model.dat))]), # lags
             c(names(model.dat)[grep("EXP_RATE_AVG", names(model.dat))])) # smooths
ice.pars <- c(c(NA, names(model.dat)[grep("ICE_LAG", names(model.dat))]), # lags
              c(NA, names(model.dat)[grep("ICE_AVG", names(model.dat))])) # smooths

combos <- expand_grid(
  df = df.pars,
  large = abund.lg.pars,
  small = abund.sm.pars,
  ice = ice.lags) %>%
  filter(!(is.na(df) & is.na(large) & is.na(small) & is.na(ice)))


# Small  ----
sm.fits <- pmap_dfr(
  combos,
  function(df, large, small, ice) {
    terms <- c(
      if (!is.na(df)) paste0("s(", df, ",k=4)") else NULL,
      if (!is.na(large)) paste0("s(", large, ",k=4)") else NULL,
      if (!is.na(small)) paste0("s(", small, ",k=4)") else NULL,
      if (!is.na(ice)) paste0("s(", ice, ",k=4)") else NULL
    )
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    fit <- gam(fml, data = sm.dat, family = betar(link = "logit"), 
               method =  "REML")
    tibble(
      df_term = df,
      large_term = large,
      small_term = small,
      ice_term = ice,
      k_terms = length(terms),
      AIC = AIC(fit),
      GCV = fit$gcv.ubre,
      edf_total = sum(fit$edf)
    )
  }
)


# Large ----
lg.fits <- pmap_dfr(
  combos,
  function(df, large, small, ice) {
    terms <- c(
      if (!is.na(df)) paste0("s(", df, ",k=4)") else NULL,
      if (!is.na(large)) paste0("s(", large, ",k=4)") else NULL,
      if (!is.na(small)) paste0("s(", small, ",k=4)") else NULL,
      if (!is.na(ice)) paste0("s(", ice, ",k=4)") else NULL
    )
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    fit <- gam(fml, data = lg.dat, family = betar(link = "logit"), 
               method =  "REML")
    tibble(
      df_term = df,
      large_term = large,
      small_term = small,
      ice_term = ice,
      k_terms = length(terms),
      AIC = AIC(fit),
      GCV = fit$gcv.ubre,
      edf_total = sum(fit$edf)
    )
  }
)




# Bind all
sm.fits %>%
  arrange(., AIC) 


sm.mod <- gam(PROP_MATURE ~ 
                s(MAT_ABUND_LG_LAG1, k = 4)+
                s(MAT_ABUND_SM_LAG1, k = 4),
                #s(EXP_RATE_AVG3, k = 4),
              family = betar(link = "logit"),
              #method =  "ML",
              data = sm.dat)

diagnose(sm.mod)

lg.fits %>%
  arrange(., AIC)

lg.mod <- gam(PROP_MATURE ~ 
                s(EXP_RATE_LAG2, k = 4)+
                s(MAT_ABUND_LG_LAG2, k = 4)+
                s(MAT_ABUND_SM_AVG2, k = 4)+
                s(ICE_LAG3, k = 4),
              family = betar(link = "logit"),
              #method =  "ML",
              data = lg.dat)

diagnose(lg.mod)

# SAM ----
sm.bioabund <-  crabpack::calc_bioabund(crab_data = ba.dat, species = "SNOW", years = years, 
                                        size_min = NULL, size_max = 95,  sex = "male", 
                                        shell_condition = "new_hardshell") %>% # all mature biomass abundance
  mutate(MAT_ABUND_SM = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000)

lg.bioabund <- crabpack::calc_bioabund(crab_data = ba.dat, species = "SNOW", years = years, 
                                       size_min = 95, size_max = NULL,  sex = "male", 
                                       shell_condition = "new_hardshell") %>% # all mature biomass abundance
  mutate(MAT_ABUND_LG = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000)

SAM.dat <-model.dat %>%
  ungroup() %>%
  dplyr::select(!contains(c("MAT_ABUND_SM", "MAT_ABUND_LG", "BIOMASS"))) %>%
  dplyr::select(!c(BIN, COUNT, MAT_COUNT, PROP_MATURE)) %>%
  distinct() %>%
  right_join(sm.bioabund %>% dplyr::select(YEAR, MAT_ABUND_SM),
             by = "YEAR") %>%
  right_join(lg.bioabund %>% dplyr::select(YEAR, MAT_ABUND_LG),
             by = "YEAR") %>%
  mutate(MAT_ABUND_SM_LAG1 = lag(MAT_ABUND_SM, n=1),
         MAT_ABUND_SM_LAG2 = lag(MAT_ABUND_SM, n=2),
         MAT_ABUND_LG_LAG1 = lag(MAT_ABUND_LG, n=1),
         MAT_ABUND_LG_LAG2 = lag(MAT_ABUND_LG, n=2),
         MAT_ABUND_SM_AVG2 = zoo::rollmean(MAT_ABUND_SM, k = 2, fill = NA, align = "center"),
         MAT_ABUND_SM_AVG3 = zoo::rollmean(MAT_ABUND_SM, k = 3, fill = NA, align = "center"),
         MAT_ABUND_LG_AVG2 = zoo::rollmean(MAT_ABUND_LG, k = 2, fill = NA, align = "center"),
         MAT_ABUND_LG_AVG3 = zoo::rollmean(MAT_ABUND_LG, k = 3, fill = NA, align = "center"))


response <- "SAM"

# SAM lags ----
# Lags and smooths
abund.lg.pars <- c(c(NA, names(SAM.dat)[grep("MAT_ABUND_LG_LAG", names(SAM.dat))]), # lags
                   c(names(SAM.dat)[grep("MAT_ABUND_LG_AVG", names(SAM.dat))])) # smooths
abund.sm.pars <- c(c(NA, names(SAM.dat)[grep("MAT_ABUND_SM_LAG", names(SAM.dat))]), # lags
                   c(names(SAM.dat)[grep("MAT_ABUND_SM_AVG", names(SAM.dat))])) # smooths
df.pars <- c(c(NA,names(SAM.dat)[grep("EXP_RATE_LAG", names(SAM.dat))]), # lags
             c(names(SAM.dat)[grep("EXP_RATE_AVG", names(SAM.dat))])) # smooths
ice.pars <- c(c(NA, names(SAM.dat)[grep("ICE_LAG", names(SAM.dat))]), # lags
              c(NA, names(SAM.dat)[grep("ICE_AVG", names(SAM.dat))])) # smooths

SAM.combos <- expand_grid(
  df = df.pars,
  large = abund.lg.pars,
  small = abund.sm.pars,
  ice = ice.lags) %>%
  filter(!(is.na(df) & is.na(large) & is.na(small) & is.na(ice)))


# SAM lags ----
SAM.fits <- pmap_dfr(
  SAM.combos,
  function(df, large, small, ice) {
    terms <- c(
      if (!is.na(df)) paste0("s(", df, ",k=4)") else NULL,
      if (!is.na(large)) paste0("s(", large, ",k=4)") else NULL,
      if (!is.na(small)) paste0("s(", small, ",k=4)") else NULL,
      if (!is.na(ice)) paste0("s(", ice, ",k=4)") else NULL
    )
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    fit <- gam(fml, data = SAM.dat, family = gaussian(), 
               method =  "REML")
    tibble(
      df_term = df,
      large_term = large,
      small_term = small,
      ice_term = ice,
      k_terms = length(terms),
      AIC = AIC(fit),
      GCV = fit$gcv.ubre,
      edf_total = sum(fit$edf)
    )
  }
)


SAM.fits %>%
  arrange(., AIC)


SAM.mod <- gam(SAM ~
                 s(EXP_RATE_AVG3, k = 4) +
                 s(MAT_ABUND_LG_LAG1, k = 4)+
                 s(MAT_ABUND_SM_AVG3, k = 4)+
                 s(ICE_LAG3, k = 4),
               family = gaussian(),
               #method =  "REML",
               data = SAM.dat)

diagnose(SAM.mod)

