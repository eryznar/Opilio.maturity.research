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
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
  reframe(SAM = weighted.mean(SIZE, weights = SAMPLING_FACTOR))

ggplot(fem.SAM, aes(LONGITUDE, LATITUDE, fill = SAM))+
  geom_point(shape = 21, stroke = NA)+
  scale_fill_viridis_c(option = "magma")+
  facet_wrap(~YEAR)+
  theme_bw()


# Weighted prop mature in 55-65
# Calculate weighted mean SAM for mature female
fem.pmat <- spec.dat.sel$specimen %>% 
  filter(SEX == 2, SIZE >=55 & SIZE <=65) %>%
  mutate(MATURE  = case_when(CLUTCH_SIZE >0 ~ 1,
                             TRUE ~ 0)) %>%
  dplyr::select(YEAR, SEX, SIZE, SAMPLING_FACTOR, MATURE, STATION_ID, LATITUDE, LONGITUDE) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
  reframe(TOT_CRAB = sum(SAMPLING_FACTOR),
          MATURE = sum(SAMPLING_FACTOR[MATURE == 1]),
          IMMATURE = TOT_CRAB - MATURE,
          PROP_MATURE = MATURE/TOT_CRAB,
          PROP_IMMATURE = IMMATURE/TOT_CRAB) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, PROP_MATURE, MATURE, IMMATURE) %>%
  rename(PMAT_5565 = PROP_MATURE) 

ggplot(fem.pmat, aes(LONGITUDE, LATITUDE, fill = PMAT_5565))+
  geom_point(shape = 21, stroke = NA)+
  scale_fill_viridis_c(direction = -1, option = "magma")+
  facet_wrap(~YEAR)+
  theme_bw()

dat <- right_join(fem.pmat, fem.SAM)

# All male large male abundance
cpue.lg <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 95, size_max = NULL,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
  reframe(LGMALE_CPUE = sum(CPUE)) 

# mature female abundance
cpue.matfem <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                                size_min = NULL, size_max = NULL,  sex = "female", 
                                                shell_condition = c("mature_female")) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
  reframe(MATFEM_CPUE = sum(CPUE)) 


# instar 1 abundance (30-50mm) (Sainte Marie?)
instar1.cpue <-  crabpack::calc_cpue(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 30, size_max = 40,  sex = "female", 
                                    shell_condition = c("new_hardshell")) %>%
  rename(INST1_CPUE = CPUE) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, INST1_CPUE) 


abund.dat <- right_join(cpue.lg, cpue.matfem, by = c("YEAR", "STATION_ID", "LATITUDE", "LONGITUDE")) %>%
  right_join(., instar1.cpue) 
unique(is.na(abund.dat))


# Bind with SAM
SAM.spatdf = right_join(dat, abund.dat)

ggplot(SAM.spatdf, aes(LONGITUDE, LATITUDE, fill = log(INST1_CPUE+100)))+
 geom_point(shape = 21, stroke = NA)+
  facet_wrap(~YEAR)+
  scale_fill_viridis_c()+
  theme_bw()

# Load Jan-April ice data
ice <- read.csv(paste0("./Maturity research/Output/spatial_ice_means_1980-", current.year, ".csv")) %>%
  #filter(name == "Mar-Apr ice") %>%
  group_by(year, latitude, longitude) %>%
  reframe(value = mean(value)) %>%
  dplyr::select(year, latitude, longitude, value) %>%
  rename(YEAR = year, LATITUDE = latitude, LONGITUDE = longitude, ICE = value)

# Load temperature occupied data
t_occ <- spec.dat.sel$haul %>%
  dplyr::select(YEAR, STATION_ID, GEAR_TEMPERATURE) %>%
  rename(FEM_TOCC = GEAR_TEMPERATURE)


# Bind all dataframes into df for modeling and plot
fem.model.dat <- right_join(SAM.spatdf,t_occ) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR)) %>%
  dplyr::select(!c(X, Y, geometry)) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1))) %>%
  arrange(YEAR) 


## ------------------------------------------------------------
## 2) Build running means (and keep in one object)
## ------------------------------------------------------------
max_lag <- 3

fem.model.dat2 <- fem.model.dat %>%
  arrange(YEAR, STATION_ID) %>%
  mutate(
    # 2‑ and 3‑year running means
    # ICE_avg2        = zoo::rollmean(ICE,         k = 2, fill = NA, align = "right"),
    # ICE_avg3        = zoo::rollmean(ICE,         k = 3, fill = NA, align = "right"),
    INST1_CPUE_avg2= zoo::rollmean(INST1_CPUE, k = 2, fill = NA, align = "right"),
    INST1_CPUE_avg3= zoo::rollmean(INST1_CPUE, k = 3, fill = NA, align = "right"),
    LGMALE_CPUE_avg2   = zoo::rollmean(LGMALE_CPUE,    k = 2, fill = NA, align = "right"),
    LGMALE_CPUE_avg3   = zoo::rollmean(LGMALE_CPUE,    k = 3, fill = NA, align = "right"),
    FEM_TOCC_avg2       = zoo::rollmean(FEM_TOCC,        k = 2, fill = NA, align = "right"),
    FEM_TOCC_avg3       = zoo::rollmean(FEM_TOCC,        k = 3, fill = NA, align = "right"),
    MATFEM_CPUE_avg2 = zoo::rollmean(MATFEM_CPUE,    k = 2, fill = NA, align = "right"),
    MATFEM_CPUE_avg3 = zoo::rollmean(MATFEM_CPUE,    k = 3, fill = NA, align = "right")
  )


# CCF ----
dat_ccf <- fem.model.dat2
vars    <- names(dat_ccf)[!names(dat_ccf) %in% c("YEAR", "YEAR_F", "SAM", "PMAT_5565", "STATION_ID", "LATITUDE", "LONGITUDE", "MATURE", "IMMATURE")]
cc_df   <- data.frame()
response <- "PMAT_5565"
for (vv in vars) {
  pp <- dat_ccf[[vv]]
  
  ok <- complete.cases(dat_ccf[[response]], pp)
  if (sum(ok) < 2) next   # only need enough paired data for CCF
  
  SAM_ok <- dat_ccf[[response]][ok]
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
  ggtitle(response)+
  theme_bw() +
  scale_x_continuous(breaks = seq(-max_lag, max_lag, 1),
                     labels = seq(-max_lag, max_lag, 1)) +
  theme(panel.grid.minor.x = element_blank())

# select top |cor| for *negative* lags (covariate leads SAM)
best_lags <- long_df %>%
  filter(lag <= 0) %>%
  group_by(short_var) %>%
  slice_max(order_by = abs(cor), n = 2, with_ties = FALSE)


# Add chosen lagged covariates (covariate precedes SAM) ----
model.dat3 <- fem.model.dat2 %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SAM, PMAT_5565, dplyr::any_of(best_lags$var)) %>%
  arrange(YEAR, STATION_ID) %>%
  mutate(LGMALE_CPUE_avg3lag3 = lag(LGMALE_CPUE_avg3, 3)) %>%
  dplyr::select(YEAR, SAM, PMAT_5565, STATION_ID, LATITUDE, LONGITUDE, # NA in PMAT is where no crab were caught at that station in that size
                FEM_TOCC, FEM_TOCC_avg2,
                INST1_CPUE_avg2, INST1_CPUE,
                LGMALE_CPUE_avg3, LGMALE_CPUE_avg3lag3,
                MATFEM_CPUE_avg2, MATFEM_CPUE_avg3) 
# CV ----
k_folds <- 5

cv_rmse <- function(fml, data, k_folds = 5) {
  data  <- dplyr::arrange(data, YEAR)
  n     <- nrow(data)
  folds <- cut(seq_len(n), breaks = k_folds, labels = FALSE)
  
  errs <- numeric(k_folds)
  
  for (k in seq_len(k_folds)) {
    test_idx  <- which(folds == k)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    train_dat <- data[train_idx, , drop = FALSE]
    test_dat  <- data[test_idx,  , drop = FALSE]
    
    fit_k <- gamm(
      fml,
      data        = train_dat,
      family      = gaussian(),
      method      = "REML",
      correlation = corAR1()  # or remove if you want IID
    )
    
    p_hat <- predict(fit_k$gam, newdata = test_dat, type = "response")
    errs[k] <- sqrt(mean((test_dat$PMAT_5565 - p_hat)^2, na.rm = TRUE))
  }
  
  mean(errs)
}

safe_cv_rmse <- function(fml, data, k_folds = 5) {
  out <- try(cv_rmse(fml, data = data, k_folds = k_folds), silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

## 4. Grid of candidate covariate terms ---------------------------
response <- "PMAT_5565"

lg.pars   <- c(NA, names(model.dat3)[grep("LGMALE_CPUE",   names(model.dat3))])
mat.pars  <- c(NA, names(model.dat3)[grep("MATFEM_CPUE",  names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("INST1_CPUE",   names(model.dat3))])
tocc.pars <- c(NA, names(model.dat3)[grep("TOCC",         names(model.dat3))])

combos <- tidyr::expand_grid(
  mat  = mat.pars,
  lg   = lg.pars,
  sm   = sm.pars,
  tocc = tocc.pars
) %>%
  dplyr::filter(!(is.na(lg) & is.na(sm) & is.na(tocc) & is.na(mat)))

safe_gamm <- purrr::safely(gamm)

fits <- purrr::pmap_dfr(
  combos,
  function(lg, sm, tocc, mat) {
    
    bio_terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ",k=4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ",k=4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
      if (!is.na(mat))  paste0("s(", mat,  ",k=4)") else NULL
    )
    
    terms <- c(
      "s(LATITUDE, LONGITUDE, bs = \"tp\", k = 30)",  # spatial interaction
      bio_terms
    )
    
    fml <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
    
    fit <- safe_gamm(
      fml,
      data        = model.dat3,
      family      = gaussian(),
      method      = "REML",
      correlation = corAR1()
    )
    
    if (!is.null(fit$error)) {
      return(tibble::tibble(
        sm_term   = sm,
        mat_term  = mat,
        male_term = lg,
        tocc_term = tocc,
        k_terms   = length(terms),
        AIC       = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        edf_total = NA_real_,
        error     = conditionMessage(fit$error)
      ))
    }
    
    cv_err <- safe_cv_rmse(fml, data = model.dat3, k_folds = k_folds)
    
    tibble::tibble(
      sm_term   = sm,
      mat_term  = mat,
      male_term = lg,
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


# Fit best model
fits %>% arrange(cv_rmse, AIC)

# fit model
mod <- gamm(
  PMAT_5565 ~ 
    s(INST1_CPUE,         k = 4) +
    s(FEM_TOCC_avg2, k = 4)+
    s(MATFEM_CPUE_avg3, k = 4)+
    s(LATITUDE, LONGITUDE, bs = "tp"),
  data        = model.dat3,
  correlation = corAR1(),
  family      = gaussian(),
)



