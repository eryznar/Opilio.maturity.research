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
  reframe(FEM_INST1_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, FEM_INST1_ABUND) 

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

M <- cor(fem.model.dat %>% dplyr::select(!c(YEAR, SAM, PMAT_5565, MATURE, IMMATURE)) %>% na.omit(), use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(M,
                   type = "upper",
                   method = "square",
                   order  = "hclust",      # cluster variables
                   addCoef.col = "black") 

femdat.long <- fem.model.dat %>%
  pivot_longer(!YEAR, names_to = "Parameter", values_to = "Value") 

ggplot(femdat.long, aes(YEAR, Value))+
  geom_line()+
  geom_point()+
  facet_wrap(~Parameter, scales = "free_y")+
  theme_bw()

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
    FEM_INST1_ABUND_avg2= zoo::rollmean(FEM_INST1_ABUND, k = 2, fill = NA, align = "right"),
    #FEM_INST1_ABUND_avg3= zoo::rollmean(FEM_INST1_ABUND, k = 3, fill = NA, align = "right"),
    MALE_LG_ABUND_avg2   = zoo::rollmean(MALE_LG_ABUND,    k = 2, fill = NA, align = "right"),
    #MALE_LG_ABUND_avg3   = zoo::rollmean(MALE_LG_ABUND,    k = 3, fill = NA, align = "right"),
    FEM_TOCC_avg2       = zoo::rollmean(FEM_TOCC,        k = 2, fill = NA, align = "right"),
    #FEM_TOCC_avg3       = zoo::rollmean(FEM_TOCC,        k = 3, fill = NA, align = "right"),
    FEM_MAT_ABUND_avg2 = zoo::rollmean(FEM_MAT_ABUND,    k = 2, fill = NA, align = "right"),
    #FEM_MAT_ABUND_avg3 = zoo::rollmean(FEM_MAT_ABUND,    k = 3, fill = NA, align = "right")
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
    
    FEM_INST1_ABUND = FEM_INST1_ABUND,
    
    FEM_MAT_ABUND  = FEM_MAT_ABUND,
    # FEM_MAT_ABUND_lag1 = lag(FEM_MAT_ABUND, 1),
    # FEM_MAT_ABUND_lag2 = lag(FEM_MAT_ABUND, 2),
    #FEM_MAT_ABUND_avg2    = FEM_MAT_ABUND_avg2,
    FEM_MAT_ABUND_lag2  = lag(FEM_MAT_ABUND_avg2, 2),
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
    
    FEM_INST1_ABUND,
    
    FEM_MAT_ABUND, FEM_MAT_ABUND_lag2,
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
sm.pars   <- c(NA, names(model.dat3)[grep("INST1_ABUND",names(model.dat3))])
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

safe_gamm <- purrr::safely(gamm)

fits <- purrr::pmap_dfr(
  combos,
  function(lg, sm, ice, tocc, mat) {
    terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ",k=4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ",k=4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k=4)") else NULL,
      if (!is.na(ice))  paste0("s(", ice,  ",k=4)") else NULL,
      if (!is.na(mat))  paste0("s(", mat,  ",k=4)") else NULL
    )
    
    ## use log(SAM) as the response
    fml <- as.formula(
      paste0("log(", response, ") ~ ", paste(terms, collapse = " + "))
    )
    
    fit <- safe_gamm(
      fml, data = model.dat3, family = gaussian(),
      method = "REML", correlation = corAR1()
    )
    
    if (!is.null(fit$error)) {
      return(tibble::tibble(
        sm_term   = sm,
        mat_term  = mat,
        male_term = lg,
        ice_term  = ice,
        tocc_term = tocc,
        k_terms   = length(terms),
        AIC       = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        R2        = NA_real_,
        edf_total = NA_real_,
        ok_resid  = FALSE,
        ok_acf    = FALSE,
        error     = conditionMessage(fit$error)
      ))
    }
    
    ## out-of-time error (on log scale)
    cv_err <- safe_cv_rmse(fml, data = model.dat3, k_folds = k_folds)
    
    ## R2 from gam component
    gam_sum <- summary(fit$result$gam)
    R2      <- gam_sum$r.sq   # or gam_sum$adj.r.sq
    
    ## residual diagnostics ----
    resid_raw <- residuals(fit$result$lme, type = "normalized")
    z_resid   <- scale(resid_raw)
    ok_resid  <- !any(abs(z_resid) > 3, na.rm = TRUE)
    
    acf_obj   <- acf(resid_raw, plot = FALSE, na.action = na.pass)
    acf_vals  <- acf_obj$acf[2:6]
    ok_acf    <- all(abs(acf_vals) < 0.4, na.rm = TRUE)
    
    tibble::tibble(
      sm_term   = sm,
      mat_term  = mat,
      male_term = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      AIC       = AIC(fit$result$lme),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = cv_err,       # RMSE on log scale
      R2        = R2,
      edf_total = sum(fit$result$gam$edf),
      ok_resid  = ok_resid,
      ok_acf    = ok_acf,
      error     = NA_character_
    )
  }
)

fits_ranked <- fits %>%
  filter(is.na(error), ok_resid, ok_acf) %>%        # diagnostics gate
  mutate(
    dAIC    = AIC - min(AIC, na.rm = TRUE),
    dRMSE   = cv_rmse - min(cv_rmse, na.rm = TRUE),
    rmse_sd = sd(cv_rmse, na.rm = TRUE)
  ) %>%
  mutate(
    keep_RMSE = dRMSE <= 2 * rmse_sd,   # main criterion
  ) %>%
  filter(keep_RMSE, dAIC <= 4) %>%     # drop only clearly bad AIC
  dplyr::select(!c(ok_resid, ok_acf, error)) %>%
  arrange(cv_rmse, AIC)

write.csv(fits_ranked, "./Maturity research/Output/SNOW_female_SAM_modelselection.csv")
read.csv("./Maturity research/Output/SNOW_female_SAM_modelselection.csv")

# fit model
mod1 <- gamm(
  log(SAM) ~ 
    #s(FEM_INST1_ABUND, k = 4)+
    s(FEM_MAT_ABUND_lag2, k = 4),
    #s(MALE_LG_ABUND, k =4)+
  #s(ICE_avg2, k = 4),
    #s(FEM_TOCC, k = 4), 
  #s(YEAR),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian(),
  method = "REML"
)

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
  xlab("Value")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_female_SAM_effectplots.png", width = 8, height = 7)

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
    
    FEM_INST1_ABUND = FEM_INST1_ABUND,
    
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
    
    ICE  = ICE,
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
    
    FEM_INST1_ABUND,
    
    FEM_MAT_ABUND_lag2, FEM_MAT_ABUND_avg2, 
    #FEM_MAT_ABUND_avg2lag2, FEM_MAT_ABUND, FEM_MAT_ABUND_lag1, FEM_MAT_ABUND_avg2lag1,
    
    MALE_LG_ABUND_avg2lag1, MALE_LG_ABUND_lag2, MALE_LG_ABUND_avg2,
   # MALE_LG_ABUND, MALE_LG_ABUND_lag1, MALE_LG_ABUND_lag2, MALE_LG_ABUND_avg2, 
    
    ICE, ICE_avg2, 
    #ICE, ICE_lag2, ICE_avg2, ICE_avg2lag1,
    
    FEM_TOCC_lag1, FEM_TOCC_avg2, 
    #FEM_TOCC_avg2lag2, FEM_TOCC ,  FEM_TOCC_lag2, FEM_TOCC_avg2lag1
  )



k_folds <- 3   # kept only for interface compatibility

cv_rmse <- function(fml,
                    data,
                    k_folds   = 3,   # ignored; kept for compatibility
                    min_train = 8,
                    gap       = 1) { # interpreted as forecast horizon h
  
  data <- data[order(data$YEAR), ]
  data <- subset(data, TOTAL_5565 > 0)  # ensure valid binomial counts
  n    <- nrow(data)
  
  h <- gap  # forecast horizon in years
  
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
    
    # fit binomial GAM on counts
    fit_k <- gam(
      fml,
      data   = train_dat,
      family = quasibinomial(link = "logit"),
      method = "REML"
    )
    
    # predict proportion matured, then compute RMSE on that scale
    p_hat <- predict(fit_k, newdata = test_dat, type = "response")
    p_obs <- with(test_dat, MATURE / TOTAL_5565)
    
    errs <- c(
      errs,
      sqrt(mean((p_obs - p_hat)^2, na.rm = TRUE))
    )
  }
  
  mean(errs, na.rm = TRUE)
}

safe_cv_rmse <- function(fml,
                         data,
                         k_folds = 3,
                         gap     = 1,
                         min_train = 8) {
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

## 4. Grid of candidate covariate terms --------------------------
response_counts <- "cbind(MATURE, TOTAL_5565 - MATURE)"

lg.pars   <- c(NA, names(model.dat3)[grep("LG_ABUND",   names(model.dat3))])
mat.pars  <- c(NA, names(model.dat3)[grep("MAT_ABUND",  names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("INST1_ABUND",names(model.dat3))])
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

k_folds <- 4
gap     <- 1

fits <- purrr::pmap_dfr(
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
        AIC       = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        edf_total = NA_real_,
        error     = conditionMessage(fit$error)
      ))
    }
    
    cv_err <- safe_cv_rmse(fml, data = model.dat3, k_folds = k_folds, gap = gap)
    
    tibble::tibble(
      sm_term   = sm,
      mat_term  = mat,
      male_term = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      AIC       = AIC(fit$result),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = cv_err,
      edf_total = sum(fit$result$edf),
      error     = NA_character_
    )
  }
)

# Evaluate model fits based on AIC and RMSE ----
fits_ranked <- fits %>%
  # remove models that failed
  dplyr::filter(is.na(error)) %>%
  # compute deltas and RMSE SD
  dplyr::mutate(
    dGCV   = GCV - min(GCV, na.rm = TRUE),
    dRMSE  = cv_rmse - min(cv_rmse, na.rm = TRUE),
    rmse_sd = sd(cv_rmse, na.rm = TRUE)
  ) %>%
  # define "neighborhood" criteria
  dplyr::mutate(
    keep_GCV  = dGCV  <= 0.7,         # e.g. within 0.01 of best GCV (tune as needed)
    keep_RMSE = dRMSE <= rmse_sd * 2   # within 2 SD of best RMSE
  ) %>%
  # keep models that are good on at least one criterion
  filter(keep_RMSE, dGCV <= 0.7) %>%    
  # rank primarily by RMSE, secondarily by GCV
  dplyr::arrange(cv_rmse, GCV)

fits_ranked

write.csv(fits_ranked, "./Maturity research/Output/SNOW_female_pmat5565_modelselection.csv")
read.csv("./Maturity research/Output/SNOW_female_pmat5565_modelselection.csv")

# fit model with temp/ice affecting recruitment ----
mod1 <- gam(
  cbind(MATURE, IMMATURE) ~ 
      s(MALE_LG_ABUND_avg2lag1, k =4) + 
    s(ICE_avg2, k = 4),
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
  facet_wrap(~ .smooth, scales = "free_x", nrow = 2) +   # facet by smooth term name
  theme_bw()+
  ylab("Partial effect")+
  xlab("Value")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_female_pmat5565_effectplots.png", width = 8, height = 7)
