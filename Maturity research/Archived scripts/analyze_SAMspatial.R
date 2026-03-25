spSAM <- read.csv("./Maturity research/Output/SAM_spatial.csv") %>%
            dplyr::select(!c(X, SAM_WTD)) %>%
            rename(SAM = SAM_UNWTD) %>%
            distinct() %>%
            filter(is.na(SAM) == FALSE) 
haul <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$haul %>%
          dplyr::select(YEAR, STATION_ID, MID_LATITUDE, MID_LONGITUDE, GEAR_TEMPERATURE, BOTTOM_DEPTH)

spSAM2 <- left_join(spSAM, haul) %>%
            filter(is.na(GEAR_TEMPERATURE) == FALSE) %>%
          dplyr::select(!c(LATITUDE, LONGITUDE, GEAR_TEMPERATURE, BOTTOM_DEPTH)) %>%
          rename(LATITUDE = MID_LATITUDE, LONGITUDE = MID_LONGITUDE)

plot.dat <- right_join(spSAM2, readRDS("./Maturity research/Output/sdmTMB_specdat.csv") %>% dplyr::select(!c(LATITUDE, LONGITUDE)), by = c("YEAR", "STATION_ID"))

# Plot
ggplot(spSAM2, aes(LONGITUDE, LATITUDE, fill = SAM))+
  geom_point(shape = 21, stroke = NA)+
  theme_bw()+
  facet_wrap(~YEAR)+
  scale_fill_gradient2(midpoint = mean(spSAM2$SAM))

# Selectivity
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  rename(SEL = selectivity, SIZE_5MM = size)

s.gam <- gam(SEL ~ s(SIZE_5MM), data = sel, family = Gamma(link = "log"))

# Specimen data accounting for selectivity
spec.mat.sel <- readRDS("./Maturity research/Output/sdmTMB_specdat.csv") %>%
  mutate(SEL = predict(s.gam, newdata= ., type = "response"), # predict size-specific selectivity
         SAMPLING_FACTOR_SEL = SAMPLING_FACTOR/SEL, # account for size specific selectivity in abundance
         SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR_SEL * PROP_MATURE,  
         SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR_SEL - SAMPLING_FACTOR_MATURE)



# Immature abundance with selectivity
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.mat.sel %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

cpue.imm.sel <-  crabpack::calc_cpue(crab_data = mat.dat.sel, species = "SNOW", 
                                                  size_min = 40, size_max = 95,  sex = "male", 
                                                  shell_condition = c("new_hardshell", "oldshell", "very_oldshell"), years = years) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(IMM_CPUE = sum(CPUE))

# All male small crab abundance
cpue.sm.sel <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 40, size_max = 95,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(SM_CPUE = sum(CPUE))


# Large male mature abundance
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.mat.sel %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

cpue.mat.sel <-  crabpack::calc_cpue(crab_data = mat.dat.sel, species = "SNOW", 
                                               size_min = 95, size_max = NULL,  sex = "male", 
                                               shell_condition = c("new_hardshell", "oldshell", "very_oldshell"), years = years) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(MAT_CPUE = sum(CPUE))

# All male large male abundance
cpue.lg.sel <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 95, size_max = NULL,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(LG_CPUE = sum(CPUE))


# All male large male abundance
cpue.male.sel <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                                  size_min = NULL, size_max = NULL,  sex = "male", 
                                                  shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(MALE_CPUE = sum(CPUE))

cpue.dat <- full_join(cpue.imm.sel, cpue.sm.sel) %>%
  full_join(., cpue.lg.sel) %>%
  full_join(., cpue.mat.sel) %>%
  full_join(., cpue.male.sel) %>%
  filter(YEAR >=1989) %>%
  mutate(PROP_SM = SM_CPUE/MALE_CPUE,
         PROP_LG = LG_CPUE/MALE_CPUE,
         PROP_MAT = MAT_CPUE/LG_CPUE) %>%
  right_join(., spec.dat.sel$haul %>% dplyr::select(YEAR, STATION_ID, GEAR_TEMPERATURE, BOTTOM_DEPTH) %>% filter(YEAR>=1989)) %>%
  full_join(., expand.grid(YEAR = 2020, STATION_ID = unique(.$STATION_ID))) %>%
  rename(BT= GEAR_TEMPERATURE) %>%
  dplyr::select(!c(PROP_SM, PROP_LG, PROP_MAT, MAT_CPUE, IMM_CPUE)) %>%
  arrange(.,  STATION_ID, YEAR) %>% # order by station within year
  group_by(STATION_ID) %>%
  mutate(SM_CPUE_AVG2 = zoo::rollmean(SM_CPUE, k = 2, fill = NA, align = "right"),
         LG_CPUE_AVG2 = zoo::rollmean(LG_CPUE, k = 2, fill = NA, align = "right"),
         BT_AVG2 = zoo::rollmean(BT, k = 2, fill = NA, align = "right"),
         SM_CPUE_LAG1 = lag(SM_CPUE, n = 1),
         LG_CPUE_LAG1 = lag(LG_CPUE, n = 1),
         BT_LAG1 = lag(BT, n = 1))
  
# Bind with SAM
SAM.cpue <- spSAM2 %>%
  right_join(cpue.dat, ., by = c("YEAR", "STATION_ID", "LATITUDE", "LONGITUDE")) # there should be no NAs in SAM here
  # at this point I don't think I need to add in missing years (really just 2020) because that is already accounted for 
  # in the station-level running averages of covariates above

# Fit models
response <- "SAM"

model.dat <- SAM.cpue

# Generate parameter combos to fit over
lg.pars <-  c(NA, names(model.dat)[grep("LG_CPUE_LAG", names(model.dat))])
sm.pars <- c(NA, names(model.dat)[grep("SM_CPUE_LAG", names(model.dat))])
bt.pars <-  c(NA,names(model.dat)[grep("BT_LAG", names(model.dat))])
yr.pars <- c(NA, "YEAR")
sp.pars <- c(NA, "STATION_ID")



combos <- expand_grid(
  #df = df.pars,
  bt= bt.pars,
  lg = lg.pars,
  sm = sm.pars,
  yr = yr.pars,
  sp = sp.pars) %>%
  filter(!(is.na(lg) & is.na(sm) & is.na(bt) & is.na(yr) & is.na(sp)))

# Fit models over parameter combinations
fits <- pmap_dfr(
  combos,
  function(bt, lg, sm, yr, sp) {
    terms <- c(
      if (!is.na(bt)) paste0("s(", bt, ",k=4)") else NULL,
      if (!is.na(lg)) paste0("s(", lg, ",k=4)") else NULL,
      if (!is.na(sm)) paste0("s(", sm, ",k=4)") else NULL,
      if (!is.na(yr)) paste0("s(", yr, ",k=4)") else NULL,
      if (!is.na(sp)) paste0("te(LONGITUDE, LATITUDE)") else NULL
    )
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    fit <- gamm(fml, data = model.dat, family = gaussian(), 
                correlation = corAR1(form = ~ YEAR | STATION_ID),
               method =  "REML")
    tibble(
      sm_term = sm,
      lg_term = lg,
      bt_term = bt,
      yr_term = yr,
      sp_term = ifelse(is.null(sp) == FALSE, "te(LONGITUDE, LATITUDE)", NA),
      k_terms = length(terms),
      AIC = AIC(fit),
      GCV = fit$gam$gcv.ubre,
      edf_total = sum(fit$gam$edf)
    )
  }
)

fits %>%
  arrange(., AIC)

mod <- gamm(SAM ~ 
              s(SM_CPUE_LAG1, k=4)+ 
              s(LG_CPUE_LAG1, k = 4)+
              s(BT_LAG1, k = 4)+
              s(YEAR, k = 4)+ 
              te(LATITUDE, LONGITUDE, k = c(10, 10)),
            random = list(STATION_ID = ~1),
           correlation = corAR1(form = ~ YEAR | STATION_ID),
           data = model.dat,
           family = gaussian())



msh <- sdmTMB::make_mesh(model.dat, c("LONGITUDE","LATITUDE"), n_knots = 100, type = "kmeans")

xtra.time <- c(2008, 2012, 2014, 2016, 2020) # missing years across all size bins

mod_dat2 <- model.dat %>%
  filter(
    is.finite(SAM),
    is.finite(SM_CPUE_LAG1),
    is.finite(LG_CPUE_LAG1),
    is.finite(BT_LAG1),
    is.finite(YEAR)
  )

msh <- sdmTMB::make_mesh(mod_dat2, c("LONGITUDE","LATITUDE"), n_knots = 100, type = "kmeans")
mod_st <- sdmTMB(
  SAM ~
    s(SM_CPUE_AVG2, k = 4) +
    s(LG_CPUE_AVG2, k = 4) +
    #s(BT_AVG2,     k = 4)+
   s(YEAR,        k = 4),
  data          = mod_dat2,
  mesh          = msh,
  spatial       = "on",
  spatiotemporal = "iid",   # IID random field over space–time
  time          = "YEAR",   # index column
  family        = gaussian(),
  anisotropy    = FALSE
)

visreg::visreg(mod_st, "SM_CPUE_AVG2")
