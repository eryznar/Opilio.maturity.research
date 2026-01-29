## ------------------------------------------------------------
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

## ------------------------------------------------------------
## 3) CCF diagnostics: which covariates lead SAM?
## ------------------------------------------------------------
dat_ccf <- model.dat2
vars    <- names(dat_ccf)[!names(dat_ccf) %in% c("YEAR", "SAM")]
cc_df   <- data.frame()

for (vv in vars) {
  pp <- dat_ccf[[vv]]
  
  ok <- complete.cases(dat_ccf$SAM, pp)
  if (sum(ok) < 2) next   # only need enough paired data for CCF
  
  SAM_ok <- dat_ccf$SAM[ok]
  pp_ok  <- pp[ok]
  
  # CCF on original (non pre-whitened) series: covariate first, response second
  cc <- ccf(pp_ok, SAM_ok, lag.max = max_lag, plot = FALSE)
  
  df <- data.frame(
    var = vv,
    lag = as.numeric(cc$lag),
    cor = as.numeric(cc$acf)
  )
  
  cc_df <- rbind(cc_df, df)
}

# Tag running‑mean variants
long_df <- cc_df %>%
  mutate(
    smooth = case_when(
      grepl("avg2", var, ignore.case = TRUE) ~ "2-year",
      grepl("avg3", var, ignore.case = TRUE) ~ "3-year",
      TRUE                                   ~ "none"
    ),
    short_var = case_when(
      grepl("avg2", var, ignore.case = TRUE) ~ gsub("_avg2", "", var, ignore.case = TRUE),
      grepl("avg3", var, ignore.case = TRUE) ~ gsub("_avg3", "", var, ignore.case = TRUE),
      TRUE                                   ~ var
    )
  )

# Plot
ggplot(long_df,
       aes(lag, cor, fill = factor(smooth, levels = c("none", "2-year", "3-year")))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("cadetblue", "salmon", "darkgoldenrod"), name = "smooth") +
  facet_wrap(~ short_var) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-max_lag, max_lag, 1),
                     labels = seq(-max_lag, max_lag, 1)) +
  theme(panel.grid.minor.x = element_blank())

# select top |cor| for *negative* lags (covariate leads SAM)
lag0 <- long_df %>%
  filter(lag == 0) %>%
  group_by(short_var) %>%
  slice_max(order_by = abs(cor), n = 1, with_ties = FALSE)

# Add one best negative lag candidate (delayed effect)
lag_neg <- long_df %>%
  filter(lag < 0) %>%
  group_by(short_var) %>%
  slice_max(order_by = abs(cor), n = 1, with_ties = FALSE)

best_lags <- bind_rows(lag0, lag_neg) %>% ungroup() %>% arrange(short_var, -abs(cor))


## ------------------------------------------------------------
## 4) Add chosen lagged covariates (covariate precedes SAM)
## ------------------------------------------------------------
model.dat3 <- model.dat2 %>%
  dplyr::select(YEAR, SAM, dplyr::any_of(best_lags$var)) %>%
  arrange(YEAR) %>%
  mutate(
    ICE_lag1 = lag(ICE, 1),
    INST1_ABUND_lag1 = lag(INST1_ABUND, 1),
    LG_ABUND_avg3lag5 = lag(LG_ABUND_avg3, 5),
    TOCC_avg3lag5 = lag(TOCC_avg3, 5)) %>%
  dplyr::select(YEAR, SAM, ICE, ICE_lag1, INST1_ABUND_avg2, INST1_ABUND_lag1, LG_ABUND, LG_ABUND_avg3lag5, 
                TOCC, TOCC_avg3lag5)

## ------------------------------------------------------------
## 5) CV function
## ------------------------------------------------------------
k_folds <- 5

cv_rmse <- function(fml, data, k_folds = 5) {
  data <- data %>% arrange(YEAR)
  n    <- nrow(data)
  folds <- cut(seq_len(n), breaks = k_folds, labels = FALSE)
  
  errs <- numeric(k_folds)
  
  for (k in seq_len(k_folds)) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    train_dat <- data[train_idx, , drop = FALSE]
    test_dat  <- data[test_idx,  , drop = FALSE]
    
    fit_k <- gamm(
      fml,
      data = train_dat,
      family = gaussian(),
      method = "REML",
      correlation = corAR1()
    )
    
    pred <- predict(fit_k$gam, newdata = test_dat, type = "response")
    errs[k] <- sqrt(mean((test_dat$SAM - pred)^2, na.rm = TRUE))
  }
  
  mean(errs)
}

## ------------------------------------------------------------
## 6) Define parameter grids on model.dat3 and run CV
## ------------------------------------------------------------
response <- "SAM"

lg.pars   <- c(NA, names(model.dat3)[grep("LG_ABUND",   names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("INST1_ABUND",names(model.dat3))])
tocc.pars <- c(NA, names(model.dat3)[grep("TOCC",       names(model.dat3))])
ice.pars  <- c(NA, names(model.dat3)[grep("ICE",        names(model.dat3))])

combos <- tidyr::expand_grid(
  ice  = ice.pars,
  lg   = lg.pars,
  sm   = sm.pars,
  tocc = tocc.pars
) %>%
  dplyr::filter(!(is.na(lg) & is.na(sm) & is.na(ice) & is.na(tocc)))

safe_gamm <- purrr::safely(gamm)

fits <- purrr::pmap_dfr(
  combos,
  function(lg, sm, ice, tocc) {
    terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ",k=4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ",k=4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
      if (!is.na(ice))  paste0("s(", ice,  ",k=4)") else NULL
    )
    
    fml <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
    
    fit <- safe_gamm(
      fml, data = model.dat3, family = gaussian(),
      method = "REML", correlation = corAR1()
    )
    
    if (!is.null(fit$error)) {
      return(tibble::tibble(
        sm_term   = sm,
        lg_term   = lg,
        ice_term  = ice,
        tocc_term = tocc,
        k_terms   = length(terms),
        AIC       = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        edf_total = NA_real_,
        error     = conditionMessage(fit$error)
      ))
    }
    
    cv_err <- cv_rmse(fml, data = model.dat3, k_folds = k_folds)
    
    tibble::tibble(
      sm_term   = sm,
      lg_term   = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      AIC       = AIC(fit$result$lme),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = cv_err,
      edf_total = sum(fit$result$gam$edf),
      error     = NA_character_
    )
  }
)

fits %>% arrange(cv_rmse, AIC)

# fit model
mod <- gamm(
  SAM ~ s(INST1_ABUND_avg2, k = 4) +
    s(LG_ABUND_avg3lag5,    k = 4)+
    s(TOCC, k = 4),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian()
)

diagnose.gamm(mod)


# Plot interaction surface
mod1 <- gamm(
  SAM ~ s(INST1_ABUND_avg2,LG_ABUND_avg3lag5, k = 4) +
    s(TOCC, k = 4),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian()
)

g1 <- mod1$gam

# grid over observed range
newdat <- expand.grid(
  INST1_ABUND_avg2   = seq(min(model.dat3$INST1_ABUND_avg2, na.rm = TRUE),
                           max(model.dat3$INST1_ABUND_avg2, na.rm = TRUE),
                           length = 50),
  LG_ABUND_avg3lag5  = seq(min(model.dat3$LG_ABUND_avg3lag5, na.rm = TRUE),
                           max(model.dat3$LG_ABUND_avg3lag5, na.rm = TRUE),
                           length = 50)
)

# hold TOCC at its mean (or another value)
newdat$TOCC <- mean(model.dat3$TOCC, na.rm = TRUE)

# predict partial effect of the interaction term
# type = "terms" returns contributions of each smooth
p <- predict(g1, newdata = newdat, type = "terms", se.fit = TRUE)

# column name for the bivariate smooth term
int_name <- "s(INST1_ABUND_avg2,LG_ABUND_avg3lag5)"

newdat$fit_int <- p$fit[, int_name]

ggplot(newdat,
       aes(INST1_ABUND_avg2, LG_ABUND_avg3lag5, fill = fit_int)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Partial effect") +
  library(dplyr)
library(ggplot2)

g1 <- mod1$gam

# grid over observed range
newdat <- expand.grid(
  INST1_ABUND_avg2   = seq(min(model.dat3$INST1_ABUND_avg2, na.rm = TRUE),
                           max(model.dat3$INST1_ABUND_avg2, na.rm = TRUE),
                           length = 50),
  LG_ABUND_avg3lag5  = seq(min(model.dat3$LG_ABUND_avg3lag5, na.rm = TRUE),
                           max(model.dat3$LG_ABUND_avg3lag5, na.rm = TRUE),
                           length = 50)
)

# hold TOCC at its mean (or another value)
newdat$TOCC <- mean(model.dat3$TOCC, na.rm = TRUE)

# predict partial effect of the interaction term
# type = "terms" returns contributions of each smooth
p <- predict(g1, newdata = newdat, type = "terms", se.fit = TRUE)

# column name for the bivariate smooth term
int_name <- "s(INST1_ABUND_avg2,LG_ABUND_avg3lag5)"

newdat$fit_int <- p$fit[, int_name]

ggplot(newdat,
       aes(INST1_ABUND_avg2, LG_ABUND_avg3lag5)) +
  geom_raster(aes(fill = fit_int)) +
  geom_contour(aes(z = fit_int), color = "black", alpha = 0.6) +
  scale_fill_viridis_c(name = "Partial effect") +
  labs(x = "INST1_ABUND_avg2",
       y = "LG_ABUND_avg3lag5") +
  theme_bw()
