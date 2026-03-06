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
  rename(SEL = selectivity, SIZE_5MM = size) %>%
  filter(year != "GAM predictions")

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

SAM.dat <- read.csv("./Maturity research/Data/SNOW_femaleSAM.csv") %>%
  dplyr::select(!X) %>%
  filter(YEAR >=1989)

ggplot(SAM.dat, aes(YEAR, SAM))+
  geom_point()+
  geom_line() +
  theme_bw()+
  geom_smooth()


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
  dplyr::select(YEAR, PROP_MATURE, MATURE, IMMATURE, TOT_CRAB) %>%
  rename(PMAT_5565 = PROP_MATURE) 

dat <- right_join(fem.pmat, SAM.dat)

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
  reframe(MALE_LG_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, MALE_LG_ABUND) 

# mature female abundance
bioabund.matfem.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = NULL, size_max = NULL,  sex = "female", 
                                            shell_condition = c("mature_female")) %>%
  group_by(YEAR) %>%
  reframe(FEM_MAT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, FEM_MAT_ABUND) 


# instar 1 abundance (30-50mm) (Sainte Marie?)
instar1 <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 35, size_max = 45,  sex = "female", 
                                    shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(FEM_COHORT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, FEM_COHORT_ABUND) 

abund.dat <- right_join(bioabund.lg.sel, bioabund.matfem.sel) %>%
  right_join(., instar1) 
unique(is.na(abund.dat))


# Bind with SAM
SAM.df = right_join(dat, abund.dat) %>%
  full_join(., data.frame(YEAR = 2020))

# Load Jan-April ice data
ice <- read.csv(paste0("./Maturity research/Output/ebs_ice_means_1980-", current.year, ".csv")) %>%
  #filter(name == "Mar-Apr ice") %>%
  group_by(year) %>%
  reframe(value = mean(value)) %>%
  dplyr::select(year, value) %>%
  rename(YEAR = year, ICE = value)

# Load temperature occupied data
t_occ <- read.csv("./Maturity research/Data/BT_occupied_females.csv") %>%
  rename(FEM_TOCC = temp_occ)

# Bind all dataframes into df for modeling and plot
fem.model.dat <- right_join(SAM.df, ice) %>%
  right_join(., t_occ %>% dplyr::select(!X)) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1))) %>%
  arrange(YEAR) %>%
  filter(YEAR >=1989) %>%
  dplyr::select(!c(SPECIES, DISTRICT))

write.csv(fem.model.dat, "./Maturity research/Output/SNOW_female_modeldata.csv")

M <- cor(fem.model.dat %>% dplyr::select(!c(YEAR, SAM, PMAT_5565, MATURE, IMMATURE)) %>% na.omit(), use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(M,
                   type = "upper",
                   method = "square",
                   order  = "hclust",      # cluster variables
                   addCoef.col = "black") 

femdat.long <- fem.model.dat %>%
  dplyr::select(!c(PMAT_5565, MATURE, IMMATURE, TOT_CRAB, SAM)) %>%
  rename("Large male abundance (≥95mm)" = "MALE_LG_ABUND",
         "Mature female abundance" = "FEM_MAT_ABUND",
         "Female cohort abundance (35-45mm)" = "FEM_COHORT_ABUND",
         "Ice % cover" = "ICE",
         "Temperature occupied" = "FEM_TOCC") %>%
  pivot_longer(!YEAR, names_to = "Parameter", values_to = "Value") 

ggplot(femdat.long, aes(YEAR, Value))+
  geom_line()+
  geom_point()+
  facet_wrap(~Parameter, scales = "free_y", ncol = 2)+
  theme_bw()+
  xlab("Year")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13))

ggsave("./Maturity research/Figures/SNOW_female_analysis_TS.png", width = 8, height = 6)


## ------------------------------------------------------------
## 2) Build running means (and keep in one object)
## ------------------------------------------------------------
max_lag <- 3

fem.model.dat2 <- fem.model.dat %>%
  arrange(YEAR) %>%
  mutate(
    # 2‑ and 3‑year running means
    ICE_avg2        = zoo::rollmean(ICE,         k = 2, fill = NA, align = "right"),
    #ICE_avg3        = zoo::rollmean(ICE,         k = 3, fill = NA, align = "right"),
    FEM_COHORT_ABUND_avg2= zoo::rollmean(FEM_COHORT_ABUND, k = 2, fill = NA, align = "right"),
    #FEM_COHORT_ABUND_avg3= zoo::rollmean(FEM_COHORT_ABUND, k = 3, fill = NA, align = "right"),
    MALE_LG_ABUND_avg2   = zoo::rollmean(MALE_LG_ABUND,    k = 2, fill = NA, align = "right"),
    #MALE_LG_ABUND_avg3   = zoo::rollmean(MALE_LG_ABUND,    k = 3, fill = NA, align = "right"),
    FEM_TOCC_avg2       = zoo::rollmean(FEM_TOCC,        k = 2, fill = NA, align = "right"),
    #FEM_TOCC_avg3       = zoo::rollmean(FEM_TOCC,        k = 3, fill = NA, align = "right"),
    FEM_MAT_ABUND_avg2 = zoo::rollmean(FEM_MAT_ABUND,    k = 2, fill = NA, align = "right"),
    #FEM_MAT_ABUND_avg3 = zoo::rollmean(FEM_MAT_ABUND,    k = 3, fill = NA, align = "right")
  )

cors <- fem.model.dat2 %>%
  dplyr::select(!c(YEAR, SAM, PMAT_5565, TOT_CRAB, MATURE, IMMATURE)) %>%
  mutate(
    MALE_LG_ABUND_lag1 = lag(MALE_LG_ABUND, 1),
    MALE_LG_ABUND_lag2 = lag(MALE_LG_ABUND, 2),
    MALE_LG_ABUND_avg2lag1 = lag(MALE_LG_ABUND_avg2, 1),
    FEM_MAT_ABUND_lag1 = lag(FEM_MAT_ABUND, 1),
    FEM_MAT_ABUND_lag2 = lag(FEM_MAT_ABUND, 2),
    FEM_MAT_ABUND_avg2lag1 = lag(FEM_MAT_ABUND_avg2, 1),
    FEM_COHORT_ABUND_lag1 = lag(FEM_COHORT_ABUND, 1),
    FEM_COHORT_ABUND_lag2 = lag(FEM_COHORT_ABUND, 2),
    FEM_COHORT_ABUND_avg2lag1 = lag(FEM_COHORT_ABUND_avg2, 1),
    ICE_lag1 = lag(ICE, 1),
    ICE_lag2 = lag(ICE, 2),
    ICE_avg2lag1 = lag(ICE_avg2, 1),
    FEM_TOCC_lag1 = lag(FEM_TOCC, 1),
    FEM_TOCC_lag2 = lag(FEM_TOCC, 2),
    FEM_TOCC_avg2lag1 = lag(FEM_TOCC_avg2, 1),
  ) %>%
  sort(colnames(.))

M <- cor(cors %>% na.omit(), use = "pairwise.complete.obs", method = "pearson") 
  
M <- M[ , sort(colnames(M))]      # reorder columns
M <- M[sort(rownames(M)), ]   
 
corrplot::corrplot(
  M,
  type      = "upper",
  method    = "color",      # or "square"
  order     = "alphabet",
  tl.col    = "black",      # label color
  tl.cex    = 0.6,          # label size
  tl.srt    = 45,           # label rotation
  addCoef.col = NA,         # no numbers on the plot
  number.cex  = 0.4,        # (used if you keep numbers)
  mar = c(0,0,1,0)          # smaller margins
)


# ------------------------------------------
## SAM 
# ------------------------------------------
response <- "SAM"

# CCF ----
dat_ccf <- fem.model.dat2
vars    <- names(dat_ccf)[!names(dat_ccf) %in% c("YEAR", "SAM", "PMAT_5565", "MATURE", "IMMATURE", "TOT_CRAB")]
cc_df   <- data.frame()
for (vv in vars) {
  pp <- dat_ccf[[vv]]
  
  ok <- complete.cases(dat_ccf[[response]], pp)
  if (sum(ok) < 2) next   # only need enough paired data for CCF
  
  SAM_ok <- dat_ccf[[response]][ok]
  pp_ok  <- pp[ok]
  
  # CCF (covariate first, response second)
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
      #grepl("avg3", var, ignore.case = TRUE) ~ "3-year",
      TRUE                                   ~ "none"
    ),
    short_var = case_when(
      grepl("avg2", var, ignore.case = TRUE) ~ gsub("_avg2", "", var, ignore.case = TRUE),
      #grepl("avg3", var, ignore.case = TRUE) ~ gsub("_avg3", "", var, ignore.case = TRUE),
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

# # select top |cor| for *negative* lags (covariate leads SAM)
best_lags <- long_df %>%
  filter(lag <= 0 & lag >=-2) %>%
  group_by(short_var) %>%
  slice_max(order_by = abs(cor), n = 3, with_ties = FALSE)


# Add chosen lagged covariates (covariate precedes SAM) ----
model.dat3 <- fem.model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    FEM_COHORT_ABUND = FEM_COHORT_ABUND,
    
    FEM_MAT_ABUND  = FEM_MAT_ABUND,
    # FEM_MAT_ABUND_lag1 = lag(FEM_MAT_ABUND, 1),
    # FEM_MAT_ABUND_lag2 = lag(FEM_MAT_ABUND, 2),
    #FEM_MAT_ABUND_avg2    = FEM_MAT_ABUND_avg2,
    FEM_MAT_ABUND_avg2lag1  = lag(FEM_MAT_ABUND_avg2, 1),
    #FEM_MAT_ABUND_avg2lag2  = lag(FEM_MAT_ABUND_avg2, 2),
    
    
    MALE_LG_ABUND  = MALE_LG_ABUND,
    # MALE_LG_ABUND_lag1 = lag(MALE_LG_ABUND, 1),
    # MALE_LG_ABUND_lag2 = lag(MALE_LG_ABUND, 2),
    MALE_LG_ABUND_avg2    = MALE_LG_ABUND_avg2,
    # MALE_LG_ABUND_avg2lag1  = lag(MALE_LG_ABUND_avg2, 1),
    # MALE_LG_ABUND_avg2lag2  = lag(MALE_LG_ABUND_avg2, 2),
    
    ICE  = ICE,
    # ICE_lag1 = lag(ICE, 1),
    # ICE_lag2 = lag(ICE, 2),
    #ICE_avg2    = ICE_avg2,
    # ICE_avg2lag1  = lag(ICE_avg2, 1),
    ICE_avg2  = ICE_avg2,
    
   FEM_TOCC  = FEM_TOCC,
    # FEM_TOCC_lag1 = lag(FEM_TOCC, 1),
    FEM_TOCC_lag2 = lag(FEM_TOCC, 2),
    # FEM_TOCC_avg2    = FEM_TOCC_avg2,
    # FEM_TOCC_avg2lag1  = lag(FEM_TOCC_avg2, 1),
    #FEM_TOCC_avg2lag2  = lag(FEM_TOCC_avg2, 2),
    
    
    TOTAL_5565              = MATURE + IMMATURE
  ) %>%
  dplyr::select(
    YEAR, SAM,
    PMAT_5565, MATURE, IMMATURE, TOTAL_5565,
    
    FEM_COHORT_ABUND,
    
    FEM_MAT_ABUND, FEM_MAT_ABUND_avg2lag1,
    #FEM_MAT_ABUND_lag1, FEM_MAT_ABUND_lag2, FEM_MAT_ABUND_avg2, FEM_MAT_ABUND_avg2lag2,
    
    MALE_LG_ABUND, MALE_LG_ABUND_avg2,
    #MALE_LG_ABUND_lag1, MALE_LG_ABUND_lag2, MALE_LG_ABUND_avg2lag2,  MALE_LG_ABUND_avg2lag1,
    
    ICE, ICE_avg2, 
    #ICE_avg2lag2, ICE_lag1, ICE_lag2, ICE_avg2lag1,
    
    FEM_TOCC_lag2, FEM_TOCC
    #FEM_TOCC_lag1, FEM_TOCC_lag2, FEM_TOCC_avg2, FEM_TOCC_avg2lag1,
  )


k_folds <- 3   # kept only for interface compatibility

cv_rmse <- function(fml,
                    data,
                    k_folds   = 3,   # ignored; kept for compatibility
                    min_train = 10,
                    gap       = 1) { # interpreted as forecast horizon h
  
  data <- data[order(data$YEAR), ]
  n    <- nrow(data)
  
  h <- gap  # forecast horizon
  
  # last origin index so that we still have h observations to assess
  last_origin <- n - h
  if (last_origin <= min_train) {
    stop("Not enough data for requested min_train / gap (h).")
  }
  
  errs <- c()
  
  # expanding (stretching) window: 1:min_train, 1:(min_train+1), ..., 1:last_origin
  for (origin in seq.int(min_train, last_origin)) {
    train_dat <- data[1:origin, , drop = FALSE]
    test_dat  <- data[(origin + 1):(origin + h), , drop = FALSE]
    
    fit_k <- gamm(
      fml,
      data        = train_dat,
      family      = gaussian(),
      method      = "REML",
      correlation = corAR1()
    )
    
    pred <- predict(fit_k$gam, newdata = test_dat, type = "response")
    
    errs <- c(
      errs,
      sqrt(mean((test_dat$SAM - pred)^2, na.rm = TRUE))
    )
  }
  
  mean(errs, na.rm = TRUE)
}

safe_cv_rmse <- function(fml, data, k_folds = 3, gap = 1, min_train = 10) {
  out <- try(
    cv_rmse(
      fml,
      data      = data,
      k_folds   = k_folds,  # ignored inside cv_rmse
      min_train = min_train,
      gap       = gap
    ),
    silent = TRUE
  )
  if (inherits(out, "try-error")) NA_real_ else out
}

# response variable name in the data
response <- "SAM"

lg.pars   <- c(NA, names(model.dat3)[grep("LG_ABUND",   names(model.dat3))])
mat.pars  <- c(NA, names(model.dat3)[grep("MAT_ABUND",  names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("COHORT_ABUND",names(model.dat3))])
tocc.pars <- c(NA, names(model.dat3)[grep("TOCC",       names(model.dat3))])
ice.pars  <- c(NA, names(model.dat3)[grep("ICE",        names(model.dat3))])

combos <- tidyr::expand_grid(
  ice  = ice.pars,
  mat  = mat.pars,
  lg   = lg.pars,
  sm   = sm.pars,
  tocc = tocc.pars
) %>%
  dplyr::filter(!(is.na(lg) & is.na(sm) & is.na(ice) & is.na(tocc) & is.na(mat)))

## ------------------------------------------------------------
## 1. Fit all models once, store AICc, no CV yet
## ------------------------------------------------------------
safe_gamm <- purrr::safely(gamm)

fits_initial <- purrr::pmap_dfr(
  combos,
  function(lg, sm, ice, tocc, mat) {
    
    terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ",k=4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ",k=4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
      if (!is.na(ice))  paste0("s(", ice,  ",k=4)") else NULL,
      if (!is.na(mat))  paste0("s(", mat,  ",k=4)") else NULL
    )
    
    fml <- as.formula(
      paste(response, "~", paste(terms, collapse = " + "))
    )
    
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
        ice_term  = ice,
        tocc_term = tocc,
        k_terms   = length(terms),
        AICc      = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        R2        = NA_real_,
        edf_total = NA_real_,
        phi       = NA_real_,
        ok_resid  = FALSE,
        ok_acf    = FALSE,
        error     = conditionMessage(fit$error)
      ))
    }
    
    ## residual diagnostics
    resid_raw <- residuals(fit$result$lme, type = "normalized")
    z_resid   <- scale(resid_raw)
    ok_resid  <- !any(abs(z_resid) > 3, na.rm = TRUE)
    
    acf_obj   <- acf(resid_raw, plot = FALSE, na.action = na.pass)
    acf_vals  <- acf_obj$acf[2:6]
    ok_acf    <- all(abs(acf_vals) < 0.4, na.rm = TRUE)
    
    gam_sum <- summary(fit$result$gam)
    
    # extract AR(1) parameter phi
    phi_val <- tryCatch(
      as.numeric(coef(fit$result$lme$modelStruct$corStruct,
                      unconstrained = FALSE)),
      error = function(e) NA_real_
    )
    
    tibble::tibble(
      sm_term   = sm,
      mat_term  = mat,
      male_term = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      AICc      = MuMIn::AICc(fit$result$lme),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = NA_real_,
      R2        = gam_sum$r.sq,
      edf_total = sum(fit$result$gam$edf),
      phi       = phi_val,
      ok_resid  = ok_resid,
      ok_acf    = ok_acf,
      error     = NA_character_
    )
  }
)

## ------------------------------------------------------------
## 2. Keep only AICc‑supported models, then run CV on them
## ------------------------------------------------------------
# Gate on diagnostics and AICc
fits_AICc <- fits_initial %>%
  filter(is.na(error), ok_resid, ok_acf) %>%
  mutate(dAICc = AICc - min(AICc, na.rm = TRUE)) %>%
  filter(dAICc <= 4)    # only models within 4 AICc units

# Now compute CV RMSE ONLY for these models
fits_AICc_cv <- fits_AICc %>%
  mutate(
    cv_rmse = purrr::pmap_dbl(
      list(male_term, sm_term, ice_term, tocc_term, mat_term),
      function(lg, sm, ice, tocc, mat) {
        
        terms <- c(
          if (!is.na(lg))   paste0("s(", lg,   ",k=4)") else NULL,
          if (!is.na(sm))   paste0("s(", sm,   ",k=4)") else NULL,
          if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
          if (!is.na(ice))  paste0("s(", ice,  ",k=4)") else NULL,
          if (!is.na(mat))  paste0("s(", mat,  ",k=4)") else NULL
        )
        
        fml <- as.formula(
          paste0("log(", response, ") ~ ", paste(terms, collapse = " + "))
        )
        
        safe_cv_rmse(fml, data = model.dat3,
                     k_folds = k_folds, gap = 1, min_train = 10)
      }
    )
  )

# ------------------------------------------------------------
## 3. Final ranking: prioritize RMSE, then AICc
## ------------------------------------------------------------
fits_ranked <- fits_AICc_cv %>%
  mutate(
    dRMSE   = cv_rmse - min(cv_rmse, na.rm = TRUE),
    rmse_sd = sd(cv_rmse, na.rm = TRUE),
    keep_RMSE = dRMSE <= 3 * rmse_sd
  ) %>%
  filter(keep_RMSE) %>%
  dplyr::select(!c(ok_resid, ok_acf, error, rmse_sd, keep_RMSE)) %>%
  arrange(cv_rmse, AICc)

write.csv(fits_ranked, "./Maturity research/Output/SNOW_female_SAM_modelselection.csv")
read.csv("./Maturity research/Output/SNOW_female_SAM_modelselection.csv")

# fit model
mod1 <- gamm(
  SAM ~ 
    s(FEM_MAT_ABUND_avg2lag1, k = 4),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian(),
  method = "REML"
)
saveRDS(mod1, "./Maturity research/Models/SNOW_femaleSAM_gamm.rda")


diagnose.gamm(mod1)

# plot facetted smooths
sm.dat <- smooth_estimates(mod1) %>%
  pivot_longer(., cols = 6:ncol(.), names_to = "resp", values_to = "value")

ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
  geom_line(color = "cadetblue", linewidth = 1.25) +
  facet_wrap(~ .smooth, scales = "free_x", nrow = 2) +   # facet by smooth term name
  theme_bw()+
  ylab("Partial effect")+
  facet_wrap(
    ~ .smooth,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(c(
      "s(FEM_MAT_ABUND_lag2)"    = "Mature female abundance (2-year lag)")))+ 
  xlab("Value")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_female_SAM_effectplots.png", width = 8, height = 5)

# ------------------------------------------
## PMAT_5565
# ------------------------------------------
## 1. CCF using PMAT_5565 as response -----------------------------
response <- "PMAT_5565"

dat_ccf <- fem.model.dat2
vars    <- names(dat_ccf)[!names(dat_ccf) %in% c("YEAR", "SAM", "PMAT_5565", "MATURE", "IMMATURE", "TOT_CRAB")]
cc_df   <- data.frame()
max_lag <- 3
for (vv in vars) {
  pp <- dat_ccf[[vv]]
  
  ok <- complete.cases(dat_ccf[[response]], pp)
  if (sum(ok) < 2) next
  
  resp_ok <- dat_ccf[[response]][ok]
  pp_ok   <- pp[ok]
  
  cc <- ccf(pp_ok, resp_ok, lag.max = max_lag, plot = FALSE)
  
  df <- data.frame(
    var = vv,
    lag = as.numeric(cc$lag),
    cor = as.numeric(cc$acf)
  )
  
  cc_df <- rbind(cc_df, df)
}

long_df <- cc_df %>%
  dplyr::mutate(
    smooth = dplyr::case_when(
      grepl("avg2", var, ignore.case = TRUE) ~ "2-year",
      TRUE                                   ~ "none"
    ),
    short_var = dplyr::case_when(
      grepl("avg2", var, ignore.case = TRUE) ~ gsub("_avg2", "", var, ignore.case = TRUE),
      TRUE                                   ~ var
    )
  )

ggplot(long_df,
       aes(lag, cor, fill = factor(smooth, levels = c("none", "2-year", "3-year")))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("cadetblue", "salmon", "darkgoldenrod"), name = "smooth") +
  facet_wrap(~ short_var) +
  theme_bw() +
  ggtitle(response) +
  scale_x_continuous(breaks = seq(-max_lag, max_lag, 1),
                     labels = seq(-max_lag, max_lag, 1)) +
  theme(panel.grid.minor.x = element_blank())

ggsave("./Maturity research/Figures/SNOW_female_pmat5565_ccf.png", width = 8, height = 7)

## 2. Add chosen lagged covariates -------------------------------
  model.dat3 <- fem.model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    FEM_COHORT_ABUND = FEM_COHORT_ABUND,
    
    #FEM_MAT_ABUND  = FEM_MAT_ABUND,
    # FEM_MAT_ABUND_lag1 = lag(FEM_MAT_ABUND, 1),
    FEM_MAT_ABUND_lag2 = lag(FEM_MAT_ABUND, 2),
    FEM_MAT_ABUND_avg2    = FEM_MAT_ABUND_avg2,
    # FEM_MAT_ABUND_avg2lag1  = lag(FEM_MAT_ABUND_avg2, 1),
    # FEM_MAT_ABUND_avg2lag2  = lag(FEM_MAT_ABUND_avg2, 2),
    
    
    # MALE_LG_ABUND  = MALE_LG_ABUND,
    # MALE_LG_ABUND_lag1 = lag(MALE_LG_ABUND, 1),
    # MALE_LG_ABUND_lag2 = lag(MALE_LG_ABUND, 2),
    #MALE_LG_ABUND_avg2    = MALE_LG_ABUND_avg2,
    MALE_LG_ABUND_avg2lag1  = lag(MALE_LG_ABUND_avg2, 1),
    MALE_LG_ABUND_avg2 = MALE_LG_ABUND_avg2,
    MALE_LG_ABUND_lag2  = lag(MALE_LG_ABUND, 2),
    
    ICE_lag1  = lag(ICE, 1),
    #ICE_lag1 = lag(ICE, 1),
    #ICE_lag2 = lag(ICE, 2),
    # ICE_avg2    = ICE_avg2,
    # ICE_avg2lag1  = lag(ICE_avg2, 1),
    ICE_avg2 = ICE_avg2,
    
    #FEM_TOCC  = FEM_TOCC,
    FEM_TOCC_lag1 = lag(FEM_TOCC, 1),
    #FEM_TOCC_lag2 = lag(FEM_TOCC, 2),
    FEM_TOCC_avg2    = FEM_TOCC_avg2,
    #FEM_TOCC_avg2lag1  = lag(FEM_TOCC_avg2, 1),
    #FEM_TOCC_avg2lag2  = lag(FEM_TOCC_avg2, 2),
    
    
    TOTAL_5565              = MATURE + IMMATURE
  ) %>%
  dplyr::select(
    YEAR, SAM,
    PMAT_5565, MATURE, IMMATURE, TOTAL_5565,
    
    FEM_COHORT_ABUND,
    
    FEM_MAT_ABUND_lag2, FEM_MAT_ABUND_avg2, 
    #FEM_MAT_ABUND_avg2lag2, FEM_MAT_ABUND, FEM_MAT_ABUND_lag1, FEM_MAT_ABUND_avg2lag1,
    
    MALE_LG_ABUND_avg2lag1, MALE_LG_ABUND_lag2, MALE_LG_ABUND_avg2,
   # MALE_LG_ABUND, MALE_LG_ABUND_lag1, MALE_LG_ABUND_lag2, MALE_LG_ABUND_avg2, 
    
    ICE_lag1, ICE_avg2, 
    #ICE, ICE_lag2, ICE_avg2, ICE_avg2lag1,
    
    FEM_TOCC_lag1, FEM_TOCC_avg2, 
    #FEM_TOCC_avg2lag2, FEM_TOCC ,  FEM_TOCC_lag2, FEM_TOCC_avg2lag1
  )


## ------------------------------------------------------------
## Time‑series CV for PMAT_5565 (binomial GAM)
## ------------------------------------------------------------
k_folds_PMAT <- 4
gap_PMAT     <- 1

cv_rmse <- function(fml,
                    data,
                    k_folds   = 3,   # ignored
                    min_train = 10,
                    gap       = 1) {
  
  data <- data[order(data$YEAR), ]
  data <- subset(data, TOTAL_5565 > 0)
  n    <- nrow(data)
  
  h <- gap
  last_origin <- n - h
  if (last_origin <= min_train) {
    stop("Not enough data for requested min_train / gap (h).")
  }
  
  errs <- c()
  eps  <- 1e-6   # small bound for proportions
  
  for (origin in seq.int(min_train, last_origin)) {
    train_dat <- data[1:origin, , drop = FALSE]
    test_dat  <- data[(origin + 1):(origin + h), , drop = FALSE]
    
    fit_k <- gam(
      fml,
      data   = train_dat,
      family = quasibinomial(link = "logit"),
      method = "REML"
    )
    
    p_hat <- predict(fit_k, newdata = test_dat, type = "response")
    
    p_obs <- with(test_dat, MATURE / TOTAL_5565)
    p_obs <- pmin(pmax(p_obs, eps), 1 - eps)
    
    errs <- c(
      errs,
      sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE))
    )
  }
  
  mean(errs, na.rm = TRUE)
}

safe_cv_rmse <- function(fml,
                         data,
                         k_folds   = 3,
                         gap       = 1,
                         min_train = 10) {
  out <- try(
    cv_rmse(
      fml,
      data      = data,
      k_folds   = k_folds,
      min_train = min_train,
      gap       = gap
    ),
    silent = TRUE
  )
  if (inherits(out, "try-error")) NA_real_ else out
}

## ------------------------------------------------------------
## Grid of candidate covariate terms
## ------------------------------------------------------------
response_counts <- "cbind(MATURE, IMMATURE)"

lg.pars   <- c(NA, names(model.dat3)[grep("LG_ABUND",   names(model.dat3))])
mat.pars  <- c(NA, names(model.dat3)[grep("MAT_ABUND",  names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("COHORT_ABUND",names(model.dat3))])
tocc.pars <- c(NA, names(model.dat3)[grep("TOCC",       names(model.dat3))])
ice.pars  <- c(NA, names(model.dat3)[grep("ICE",        names(model.dat3))])

combos <- tidyr::expand_grid(
  ice  = ice.pars,
  mat  = mat.pars,
  lg   = lg.pars,
  sm   = sm.pars,
  tocc = tocc.pars
) %>%
  dplyr::filter(!(is.na(lg) & is.na(sm) & is.na(ice) & is.na(tocc) & is.na(mat)))

safe_gam <- purrr::safely(gam)

## ------------------------------------------------------------
## 1. Initial fits: store GCV, no CV yet
## ------------------------------------------------------------
fits_initial <- purrr::pmap_dfr(
  combos,
  function(lg, sm, ice, tocc, mat) {
    
    terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ", k = 4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ", k = 4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ", k = 4)") else NULL,
      if (!is.na(ice))  paste0("s(", ice,  ", k = 4)") else NULL,
      if (!is.na(mat))  paste0("s(", mat,  ", k = 4)") else NULL
    )
    
    rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
    fml <- as.formula(paste(response_counts, "~", rhs))
    
    fit <- safe_gam(
      fml,
      data   = subset(model.dat3, TOTAL_5565 > 0),
      family = quasibinomial(link = "logit"),
      method = "REML"
    )
    
    if (!is.null(fit$error)) {
      return(tibble::tibble(
        sm_term   = sm,
        mat_term  = mat,
        male_term = lg,
        ice_term  = ice,
        tocc_term = tocc,
        k_terms   = length(terms),
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        edf_total = NA_real_,
        error     = conditionMessage(fit$error)
      ))
    }
    
    tibble::tibble(
      sm_term   = sm,
      mat_term  = mat,
      male_term = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = NA_real_,
      edf_total = sum(fit$result$edf),
      error     = NA_character_
    )
  }
)

## ------------------------------------------------------------
## 2. Keep only GCV‑supported models, then run CV
## ------------------------------------------------------------
fits_GCV <- fits_initial %>%
  dplyr::filter(is.na(error)) %>%
  dplyr::mutate(
    dGCV = GCV - min(GCV, na.rm = TRUE)
  ) %>%
  dplyr::filter(dGCV <= 0.7)   # choose a small GCV window you like

# helper that returns list(RMSE = ..., converged = TRUE/FALSE)
cv_rmse_flag <- function(fml, data, k_folds, gap, min_train) {
  out <- try(
    cv_rmse(
      fml,
      data      = data,
      k_folds   = k_folds,
      min_train = min_train,
      gap       = gap
    ),
    silent = TRUE
  )
  
  if (inherits(out, "try-error") || any(grepl("did not converge", out))) {
    list(RMSE = NA_real_, converged = FALSE)
  } else {
    list(RMSE = out, converged = TRUE)
  }
}




fits_GCV_cv <- fits_GCV %>%
  mutate(
    cv_res = purrr::pmap(
      list(male_term, sm_term, ice_term, tocc_term, mat_term),
      function(lg, sm, ice, tocc, mat) {
        
        terms <- c(
          if (!is.na(lg))   paste0("s(", lg,   ", k = 4)") else NULL,
          if (!is.na(sm))   paste0("s(", sm,   ", k = 4)") else NULL,
          if (!is.na(tocc)) paste0("s(", tocc, ", k = 4)") else NULL,
          if (!is.na(ice))  paste0("s(", ice,  ", k = 4)") else NULL,
          if (!is.na(mat))  paste0("s(", mat,  ", k = 4)") else NULL
        )
        
        rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
        fml <- as.formula(paste(response_counts, "~", rhs))
        
        cv_rmse_flag(
          fml,
          data      = model.dat3,
          k_folds   = k_folds_PMAT,
          gap       = gap_PMAT,
          min_train = 10
        )
      }
    ),
    cv_rmse    = purrr::map_dbl(cv_res, "RMSE"),
    cv_ok_conv = purrr::map_lgl(cv_res, "converged")
  ) %>%
  dplyr::select(-cv_res)

## ------------------------------------------------------------
## 3. Final ranking: prioritize RMSE, then GCV
## ------------------------------------------------------------
fits_ranked <- fits_GCV_cv %>%
  dplyr::mutate(
    dRMSE   = cv_rmse - min(cv_rmse, na.rm = TRUE),
    rmse_sd = sd(cv_rmse, na.rm = TRUE),
    keep_RMSE = dRMSE <= 2 * rmse_sd
  ) %>%
  dplyr::filter(keep_RMSE) %>%
  dplyr::select(!c(error, rmse_sd, keep_RMSE)) %>%
  dplyr::arrange(cv_rmse, GCV)

write.csv(fits_ranked, "./Maturity research/Output/snow_female_pmat5565_modelselection.csv")

# Fit best model
mod1 <- gam(
  cbind(MATURE, IMMATURE) ~ 
    #s(FEM_MAT_ABUND_lag2, k =4) +
      s(MALE_LG_ABUND_avg2lag1, k =4) + 
    s(ICE_avg2, k = 4),
  data   = model.dat3,
  family = quasibinomial(link = "logit"),
  method = "REML"
)

saveRDS(mod1, "./Maturity research/Models/SNOW_femalepmat5565_gam.rda")

# Fit best model
mod1 <- gam(
  cbind(MATURE, IMMATURE) ~ 
    #s(FEM_MAT_ABUND_lag2, k =4) +
    s(MALE_LG_ABUND_avg2lag1, k =4) + 
    s(ICE_avg2, k = 4)+
    s(FEM_TOCC_avg2, k=4),
  data   = model.dat3,
  family = quasibinomial(link = "logit"),
  method = "REML"
)
gam.check(mod1)
summary(mod1)
acf(mod1$residuals)

# plot facetted smooths
sm.dat <- smooth_estimates(mod1) %>%
  pivot_longer(., cols = 6:ncol(.), names_to = "resp", values_to = "value")

ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
  geom_line(color = "cadetblue", linewidth = 1.25) +
  facet_wrap(
    ~ .smooth,
    scales = "free_x",
    nrow = 2,
    labeller = as_labeller(c(
      "s(ICE_avg2)"    = "Ice % cover (2-year avg)",
      "s(MALE_LG_ABUND_avg2lag1)"   = "Large male abundance (2-year avg, 1-year lag)"))
    )+ 
  theme_bw()+
  ylab("Partial effect")+
  xlab("Value")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_female_pmat5565_effectplots.png", width = 8, height = 7)

# Data actually used in the model (after any na.omit etc.)
dat_use <- mod1$model

# Deviance residuals in the same order/rows
res <- residuals(mod1, type = "deviance")

par(mfrow = c(1, 2))

with(dat_use, {
  plot(MALE_LG_ABUND_avg2lag1, res,
       xlab = "MALE_LG_ABUND_avg2lag1",
       ylab = "Deviance residuals")
  abline(h = 0, col = "red", lty = 2)
  
  plot(ICE_avg2, res,
       xlab = "ICE_avg2",
       ylab = "Deviance residuals")
  abline(h = 0, col = "red", lty = 2)
})

par(mfrow = c(1, 1))