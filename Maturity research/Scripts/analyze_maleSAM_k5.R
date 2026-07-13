## ------------------------------------------------------------
# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
# sdmTMB model
mod <- readRDS("./Maturity research/Models/snowmale_sdmTMB_spVAR_noBIN_k300.rda")

# SAM
SAM.dat <- read.csv("./Maturity research/Data/SNOW_maleSAM_k5.csv") %>%
  dplyr::select(!X) %>%
  full_join(., expand.grid(YEAR = 2020))


ggplot(SAM.dat, aes(YEAR, SAM))+
  geom_line()+
  geom_point()+
  theme_bw()+
  annotate(geom = "text", x=2000, y = 80, label = "p'<0.05", size = 5)+
  geom_smooth(method = "lm")

summary(lme(SAM ~ YEAR, data = na.omit(SAM.dat), random = ~ 1 | YEAR, correlation = corAR1()))


# Selectivity
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  rename(SEL = selectivity, SIZE_5MM = size) %>%
  filter(year != "GAM predictions")

s.gam <- gam(SEL ~ s(SIZE_5MM, k = 5), method = "REML", data = sel, family = Gamma(link = "log"))


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

# All male large male abundance
bioabund.lg.sel <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                            size_min = 95, size_max = NULL,  sex = "male", 
                                            shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(LG_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, LG_ABUND) %>%
  filter(YEAR >=1989)

# instar 1 abundance (40-60mm) (Sainte Marie?)
instar1 <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 40, size_max = 60,  sex = "male", 
                                    shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(COHORT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, COHORT_ABUND)  %>%
  filter(YEAR >= 1989)

abund.dat <- right_join(bioabund.lg.sel, instar1)
unique(is.na(abund.dat))


# Bind with SAM
SAM.abund <- right_join(SAM.dat, abund.dat) 

unique(is.na(SAM.abund)) # SAM should have NAs in years where no chela were measured
unique(SAM.abund[is.na(SAM.abund$SAM) == TRUE,]$YEAR) #2008, 2012, 2014, 2016

# Load Jan-April ice data
ice <- read.csv(paste0("./Maturity research/Output/ebs_ice_means_1980-", current.year, ".csv")) %>%
  #filter(name == "Mar-Apr ice") %>%
  # group_by(year) %>%
  # reframe(value = mean(value)) %>%
  dplyr::select(year, value) %>%
  rename(YEAR = year, ICE = value)

# Load temperature occupied data
t_occ <- read.csv("./Maturity research/Data/BT_occupied_k5.csv") %>%
  rename(TOCC = temp_occ)

# Bind all dataframes into df for modeling and plot
model.dat <- right_join(SAM.abund, ice) %>%
  right_join(., t_occ) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1))) %>%
  arrange(YEAR) %>%
  dplyr::select(!c(X.1,X, SPECIES, DISTRICT))



## ------------------------------------------------------------
## 2) Build running means (and keep in one object)
## ------------------------------------------------------------
max_lag <- 3

model.dat2 <- model.dat %>%
  arrange(YEAR) %>%
  mutate(
    # 2‑ and 3‑year running means
    ICE_avg2        = zoo::rollmean(ICE,         k = 2, fill = NA, align = "right"),
    #ICE_avg3        = zoo::rollmean(ICE,         k = 3, fill = NA, align = "right"),
    COHORT_ABUND_avg2= zoo::rollmean(COHORT_ABUND, k = 2, fill = NA, align = "right"),
    #COHORT_ABUND_avg3= zoo::rollmean(COHORT_ABUND, k = 3, fill = NA, align = "right"),
    LG_ABUND_avg2   = zoo::rollmean(LG_ABUND,    k = 2, fill = NA, align = "right"),
    #LG_ABUND_avg3   = zoo::rollmean(LG_ABUND,    k = 3, fill = NA, align = "right"),
    TOCC_avg2       = zoo::rollmean(TOCC,        k = 2, fill = NA, align = "right"),
    #TOCC_avg3       = zoo::rollmean(TOCC,        k = 3, fill = NA, align = "right")
  )


cors <- model.dat2 %>%
  dplyr::select(!c(YEAR, SAM))%>%
  mutate(ICE_lag1 = lag(ICE, 1),
         ICE_lag2 = lag(ICE, 2),
         ICE_avg2lag1 = lag(ICE_avg2, 1),
         COHORT_ABUND_lag1 = lag(COHORT_ABUND, 1),
         COHORT_ABUND_lag2 = lag(COHORT_ABUND, 2),
         COHORT_ABUND_avg2lag1 = lag(COHORT_ABUND_avg2, 1),
         LG_ABUND_lag1 = lag(LG_ABUND, 1),
         LG_ABUND_lag2 = lag(LG_ABUND, 2),
         LG_ABUND_avg2lag1 = lag(LG_ABUND_avg2, 1),
         TOCC_lag1 = lag(TOCC, 1),
         TOCC_lag2 = lag(TOCC, 2),
         TOCC_avg2lag1 = lag(TOCC_avg2, 1))

M <- cor(cors, use = "pairwise.complete.obs", method = "pearson")
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
ggplot(long_df %>% filter(lag<=0, !(lag == -3 & smooth == "2-year")),
       aes(lag, cor, fill = factor(smooth, levels = c("none", "2-year", "3-year")))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("cadetblue", "salmon", "darkgoldenrod"), name = "smooth") +
  facet_wrap(~ short_var) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-max_lag, max_lag, 1),
                     labels = seq(-max_lag, max_lag, 1)) +
  theme(panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size =14),
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14))

#ggsave("./Maturity research/Figures/SNOW_male_SAM_ccf.png", width = 5, height = 5)

long_df %>%
  group_by(short_var) %>%
  filter(lag <=0) %>%
  mutate(abs_cor = abs(cor)) %>%
  slice_max(order_by = abs_cor, n = 2, with_ties = FALSE)


long_df %>%
  group_by(short_var) %>%
  filter(lag <=0, lag>-3) %>%
  mutate(abs_cor = abs(cor)) %>%
  slice_min(order_by = abs_cor, n = 1, with_ties = FALSE)

## ------------------------------------------------------------
## 4) Add chosen lagged covariates (covariate precedes SAM)
## ------------------------------------------------------------
model.dat3 <- model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    COHORT_ABUND = COHORT_ABUND, # hard coding to lag0, no avg
    
    LG_ABUND  = LG_ABUND,
    #LG_ABUND_lag1 = lag(LG_ABUND, 1),
    #LG_ABUND_lag2 = lag(LG_ABUND, 2),
    LG_ABUND_avg2    = LG_ABUND_avg2,
    #LG_ABUND_avg2lag1  = lag(LG_ABUND_avg2, 1),
    #LG_ABUND_avg2lag2  = lag(LG_ABUND_avg2, 2),
    
    ICE  = ICE,
    #ICE_lag1 = lag(ICE, 1),
    #ICE_lag2 = lag(ICE, 2),
    ICE_avg2    = ICE_avg2,
    #ICE_avg2lag1  = lag(ICE_avg2, 1),
    #ICE_avg2lag2  = lag(ICE_avg2, 2),
    
    TOCC  = TOCC,
    TOCC_avg2 = TOCC_avg2,
    #TOCC_lag1 = lag(TOCC, 1),
    #TOCC_lag2 = lag(TOCC, 2),
    #TOCC_avg2lag1    = lag(TOCC_avg2, 1)
    #TOCC_avg2lag1  = lag(TOCC_avg2, 1),
    #TOCC_avg2lag2  = lag(TOCC_avg2, 2),
    
  ) %>%
  dplyr::select(
    YEAR, SAM,
    
    COHORT_ABUND,
    
    LG_ABUND, LG_ABUND_avg2,
    #LG_ABUND_avg2lag2, 
    #LG_ABUND_lag1, LG_ABUND_lag2, 
    
    ICE, ICE_avg2,
    #ICE_lag1, ICE_lag2, ICE_avg2lag2,
    
    TOCC_avg2, TOCC,
    #TOCC_lag1, TOCC_lag2, TOCC_avg2lag2
  )


worst.lags <-  model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    COHORT_ABUND_lag2 = lag(COHORT_ABUND, 2), # hard coding to lag0, no avg
    
    LG_ABUND_lag2  = lag(LG_ABUND, 2),
    
    ICE_avg2lag2  = lag(ICE_avg2, 2),
    
    TOCC_avg2lag2 = lag(TOCC_avg2, 2),
    
  ) %>%
  dplyr::select(
    YEAR, SAM,
    
    COHORT_ABUND_lag2,
    
    LG_ABUND_lag2,
    
    ICE_avg2lag2,
    TOCC_avg2lag2)


## ------------------------------------------------------------
## 5) CV function
## ------------------------------------------------------------
k_folds <- 5

cv_rmse <- function(fml,
                    data,
                    k_folds   = 5,   # kept for interface compatibility
                    min_train = 10,
                    gap       = 1,   # forecast horizon in years
                    ...) {
  
  data <- data[order(data$YEAR), ]
  n    <- nrow(data)
  
  h <- gap
  last_origin <- n - h
  if (last_origin <= min_train) {
    stop("Not enough data for requested min_train / gap (h).")
  }
  
  errs <- c()
  
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

safe_cv_rmse <- function(fml,
                         data,
                         k_folds   = 5,
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
## 6) Parameter grid and 2‑stage selection (AICc then CV)
## ------------------------------------------------------------
response <- "SAM"

lg.pars   <- c(NA, names(model.dat3)[grep("LG_ABUND",    names(model.dat3))])
sm.pars   <- c(NA, names(model.dat3)[grep("COHORT_ABUND", names(model.dat3))])
tocc.pars <- c(NA, names(model.dat3)[grep("TOCC",        names(model.dat3))])
ice.pars  <- c(NA, names(model.dat3)[grep("ICE",         names(model.dat3))])

combos <- tidyr::expand_grid(
  ice  = ice.pars,
  lg   = lg.pars,
  sm   = sm.pars,
  tocc = tocc.pars
) %>%
  dplyr::filter(!(is.na(lg) & is.na(sm) & is.na(ice) & is.na(tocc)))

safe_gamm <- purrr::safely(gamm)

fits_initial <- purrr::pmap_dfr(
  combos,
  function(lg, sm, ice, tocc) {
    
    terms <- c(
      if (!is.na(lg))   paste0("s(", lg,   ",k = 4)") else NULL,
      if (!is.na(sm))   paste0("s(", sm,   ",k = 4)") else NULL,
      if (!is.na(tocc)) paste0("s(", tocc, ",k = 4)") else NULL,
      if (!is.na(ice))  paste0("s(", ice,  ",k = 4)") else NULL
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
        lg_term   = lg,
        ice_term  = ice,
        tocc_term = tocc,
        k_terms   = length(terms),
        AICc      = NA_real_,
        GCV       = NA_real_,
        cv_rmse   = NA_real_,
        edf_total = NA_real_,
        phi       = NA_real_,
        error     = conditionMessage(fit$error)
      ))
    }
    
    # extract AR(1) parameter (phi)
    phi_val <- tryCatch(
      as.numeric(coef(fit$result$lme$modelStruct$corStruct, unconstrained = FALSE)),
      error = function(e) NA_real_
    )
    
    tibble::tibble(
      sm_term   = sm,
      lg_term   = lg,
      ice_term  = ice,
      tocc_term = tocc,
      k_terms   = length(terms),
      AICc      = MuMIn::AICc(fit$result$lme),
      GCV       = fit$result$gcv.ubre,
      cv_rmse   = NA_real_,
      edf_total = sum(fit$result$gam$edf),
      phi       = phi_val,
      error     = NA_character_
    )
  }
)

## 6b. Keep only AICc‑supported models, then run CV ----------------------
fits_AICc <- fits_initial %>%
  filter(is.na(error)) %>%
  mutate(dAICc = AICc - min(AICc, na.rm = TRUE)) %>%
  filter(dAICc <= 4)

fits_AICc_cv <- fits_AICc %>%
  mutate(
    cv_rmse = purrr::pmap_dbl(
      list(lg_term, sm_term, ice_term, tocc_term),
      function(lg, sm, ice, tocc) {
        
        terms <- c(
          if (!is.na(lg))   paste0("s(", lg,   ",k = 4)") else NULL,
          if (!is.na(sm))   paste0("s(", sm,   ",k = 4)") else NULL,
          if (!is.na(tocc)) paste0("s(", tocc, ",k = 4)") else NULL,
          if (!is.na(ice))  paste0("s(", ice,  ",k = 4)") else NULL
        )
        
        fml <- as.formula(paste(response, "~", paste(terms, collapse = " + ")))
        
        safe_cv_rmse(
          fml,
          data      = model.dat3,
          k_folds   = k_folds,
          gap       = 1,
          min_train = 10
        )
      }
    )
  )

## 6c. Final ranking: prioritize RMSE, then AICc --------------------------
fits_ranked <- fits_AICc_cv %>%
  mutate(
    dRMSE   = cv_rmse - min(cv_rmse, na.rm = TRUE),
    rmse_sd = sd(cv_rmse, na.rm = TRUE),
    keep_RMSE = dRMSE <= 2 * rmse_sd
  ) %>%
  filter(keep_RMSE) %>%
  dplyr::select(!c(error, rmse_sd, keep_RMSE)) %>%
  arrange(cv_rmse, AICc)

fits_ranked

#write.csv(fits_ranked, "./Maturity research/Output/SNOW_male_SAM_modelselection.csv")

#read.csv("./Maturity research/Output/SNOW_male_SAM_modelselection.csv")

# fit model
mod <- gamm(
  SAM ~ s(COHORT_ABUND, k = 4) +
    s(LG_ABUND_avg2,    k = 4)+
  #s(ICE_avg2lag1, k = 4),
  s(TOCC, k = 4),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian()
)

mod2 <- gamm(
  SAM ~ s(COHORT_ABUND, k = 4) +
    s(LG_ABUND_avg2,    k = 4),
    #s(ICE_avg2lag1, k = 4),
    #s(TOCC, k = 4),
  correlation = corAR1(),
  data        = model.dat3,
  family      = gaussian()
)

saveRDS(mod2, "./Maturity research/Models/SNOW_maleSAM_gamm_k5.rda")

diagnose.gamm(mod)
diagnose.gamm(mod2)

mod2 <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm_k5.rda")
# plot facetted smooths
sm.dat <- smooth_estimates(mod2) %>%
  pivot_longer(., cols = 6:ncol(.), names_to = "resp", values_to = "value") 

dat_pr <- gratia::add_partial_residuals(
  data  = mod2$gam$model,
  model = mod2$gam
)

smooth_cols <- grep("^s\\(", names(dat_pr), value = TRUE)

dat_pr_long <- lapply(smooth_cols, function(sc) {
  covar <- str_match(sc, "s\\(([^)]+)\\)")[, 2]
  tibble(
    .smooth  = sc,
    x        = dat_pr[[covar]],
    .partial = dat_pr[[sc]]
  )
}) %>%
  bind_rows()

ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
  geom_line(color = "cadetblue", linewidth = 1.25) +
  theme_bw()+
  ylab("Partial effect")+
  facet_wrap(
    ~ .smooth,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(c(
      "s(COHORT_ABUND)"    = "Male cohort size",
      "s(LG_ABUND_avg2)" = "Large male abundance (t, t-1)")))+ 
  xlab("Abundance (millions)")+
  geom_rug(data = dat_pr_long,
           aes(x = x),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.4) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        aspect.ratio = 3/4,
        plot.margin = margin(t = 5, r = 5, b = 8, l = 5, unit = "mm"))


ggsave("./Maturity research/Figures/SNOW_male_SAM_effectplots_k5.png", width = 8, height = 4)


# Fit worst SAM model ------
worst.lags

worst.mod <- gamm(
  SAM ~ s(COHORT_ABUND_lag2, k = 4) +
    s(LG_ABUND_lag2,    k = 4),
  correlation = corAR1(),
  data        = worst.lags,
  family      = gaussian()
)

saveRDS(worst.mod, "./Maturity research/Models/SNOW_maleSAM_worstmod_k5.rda")



# PROP_INDUSTRY PREFERRED ----
mod <- readRDS("./Maturity research/Models/snowmale_sdmTMB_spVAR_noBIN_k300.rda")
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")

# #Filter predicted specimen data by params (not size yet for full join)
spec.dat.mat <- spec.dat$specimen %>%
  filter(YEAR %in% mod$data$YEAR, SHELL_CONDITION == 2, SEX == 1, !c(YEAR == 2025 & SIZE == 175.9)) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2,
         YEAR_SCALED = scale(YEAR)) %>%
  mutate(SEL = predict(s.gam, newdata= ., type = "response"), # predict size-specific selectivity
         SAMPLING_FACTOR_SEL = SAMPLING_FACTOR/SEL) %>% # account for size specific selectivity in abundance
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR)) %>%
  predict(mod, ., type = "response", se = FALSE) %>%
  rename(PROP_MATURE = est) %>%
  mutate(SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR_SEL * PROP_MATURE,
         SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR_SEL-SAMPLING_FACTOR_MATURE)


# Mature abundance >=101 SH2
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.dat.mat %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

ind.pref <-  crabpack::calc_bioabund(crab_data = mat.dat.sel, species = "SNOW", 
                                     size_min = 101, size_max = NULL,  sex = "male", 
                                     shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(IND_PREF = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, IND_PREF) 

# All mature abundance SH2
mat.dat.sel <- spec.dat
mat.dat.sel$specimen <- spec.dat.mat %>%
  dplyr::select(!SAMPLING_FACTOR) %>% # removing original SF
  rename(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # renaming mature SF to SF so crabpack recognizes, this accounts for sel

all.mat <-  crabpack::calc_bioabund(crab_data = mat.dat.sel, species = "SNOW", 
                                    size_min = NULL, size_max = NULL,  sex = "male", 
                                    shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(ALL_MAT = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ALL_MAT) 

all.mat$PROP_INDPREF <- ind.pref$IND_PREF/all.mat$ALL_MAT

all.mat %>% 
  filter(YEAR >= 1989) %>%
  full_join(., data.frame(YEAR = 2020)) -> all.mat2

ggplot(all.mat2, aes(YEAR, PROP_INDPREF))+
  geom_line()+
  geom_point()+
  theme_bw()

#ggsave("./Maturity research/Figures/SNOW_male_propindpref.png", width = 9, height = 7)


indpref.dat <- right_join(ind.pref, all.mat) %>%
  filter(YEAR >= 1989) %>%
  right_join(model.dat) %>%
  full_join(expand.grid(YEAR = 2020))

write.csv(indpref.dat, "./Maturity research/Output/SNOW_male_modeldata_k5.csv")


M <- cor(indpref.dat %>% dplyr::select(!c(SAM, YEAR, IND_PREF, ALL_MAT, PROP_INDPREF)), use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(M,
                   type = "upper",
                   method = "square",
                   order  = "hclust",      # cluster variables
                   addCoef.col = "black") 


# Evaluate lags ----
max_lag <- 3

model.dat2 <- indpref.dat %>%
  dplyr::select(!SAM) %>%
  arrange(YEAR) %>%
  mutate(
    # 2‑ and 3‑year running means
    ICE_avg2        = zoo::rollmean(ICE,         k = 2, fill = NA, align = "right"),
    #ICE_avg3        = zoo::rollmean(ICE,         k = 3, fill = NA, align = "right"),
    COHORT_ABUND_avg2= zoo::rollmean(COHORT_ABUND, k = 2, fill = NA, align = "right"),
    #COHORT_ABUND_avg3= zoo::rollmean(COHORT_ABUND, k = 3, fill = NA, align = "right"),
    LG_ABUND_avg2   = zoo::rollmean(LG_ABUND,    k = 2, fill = NA, align = "right"),
    #LG_ABUND_avg3   = zoo::rollmean(LG_ABUND,    k = 3, fill = NA, align = "right"),
    TOCC_avg2       = zoo::rollmean(TOCC,        k = 2, fill = NA, align = "right"),
    #TOCC_avg3       = zoo::rollmean(TOCC,        k = 3, fill = NA, align = "right")
  )

M <- cor(model.dat2 %>% dplyr::select(!c(YEAR, IND_PREF, ALL_MAT, PROP_INDPREF)) %>% na.omit(), use = "pairwise.complete.obs", method = "pearson")
corrplot::corrplot(M,
                   type = "upper",
                   method = "square",
                   order  = "hclust",      # cluster variables
                   addCoef.col = "black") 


## ------------------------------------------------------------
## 3) CCF diagnostics
## ------------------------------------------------------------
dat_ccf <- model.dat2
vars    <- names(dat_ccf)[!names(dat_ccf) %in% c("YEAR", "PROP_INDPREF", "ALL_MAT", "IND_PREF")]
cc_df   <- data.frame()

for (vv in vars) {
  pp <- dat_ccf[[vv]]
  
  ok <- complete.cases(dat_ccf$PROP_INDPREF, pp)
  if (sum(ok) < 2) next   # only need enough paired data for CCF
  
  pind_ok <- dat_ccf$PROP_INDPREF[ok]
  pp_ok  <- pp[ok]
  
  # CCF on original (non pre-whitened) series: covariate first, response second
  cc <- ccf(pp_ok, pind_ok, lag.max = max_lag, plot = FALSE)
  
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

#ggsave("./Maturity research/Figures/SNOW_male_indpref_ccf.png", width = 8, height = 7)

long_df %>%
  group_by(short_var) %>%
  filter(lag <=0, !(lag == -3 & smooth == "2-year")) %>%
  mutate(abs_cor = abs(cor)) %>%
  slice_max(order_by = abs_cor, n = 2, with_ties = FALSE)

long_df %>%
  group_by(short_var) %>%
  filter(lag <=0, !(lag == -3 & smooth == "2-year")) %>%
  mutate(abs_cor = abs(cor)) %>%
  slice_min(order_by = abs_cor, n = 1, with_ties = FALSE)

## ------------------------------------------------------------
## 4) Add chosen lagged covariates 
## ------------------------------------------------------------
model.dat3 <- model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    COHORT_ABUND = COHORT_ABUND, # hard coding to lag0, no avg
    
    LG_ABUND_lag1  = lag(LG_ABUND, 1),
    #LG_ABUND_lag1 = lag(LG_ABUND, 1),
    #LG_ABUND_lag2 = lag(LG_ABUND, 2),
    LG_ABUND_avg2    = LG_ABUND_avg2,
    #LG_ABUND_avg2lag1  = lag(LG_ABUND_avg2, 1),
    #LG_ABUND_avg2lag2  = lag(LG_ABUND_avg2, 2),
    
    ICE_avg2  = ICE_avg2,
    #ICE_lag1 = lag(ICE, 1),
    #ICE_lag2 = lag(ICE, 2),
    ICE   = ICE,
    #ICE_avg2lag1  = lag(ICE_avg2, 1),
    #ICE_avg2lag2  = lag(ICE_avg2, 2),
    
    #TOCC_avg2  = TOCC_avg2,
    #TOCC_lag1 = lag(TOCC, 1),
    #TOCC_lag2 = lag(TOCC, 2),
    TOCC_lag1    = lag(TOCC, 1),
    TOCC_avg2lag1  = lag(TOCC_avg2, 1),
    #TOCC_avg2lag2  = lag(TOCC_avg2, 2),
    
  ) %>%
  dplyr::select(
    YEAR, PROP_INDPREF, ALL_MAT, IND_PREF,
    
    COHORT_ABUND,
    
    LG_ABUND_lag1, LG_ABUND_avg2,
    #LG_ABUND_avg2lag2, 
    #LG_ABUND_lag1, LG_ABUND_lag2, 
    
    ICE_avg2, ICE,
    #ICE_lag1, ICE_lag2, ICE_avg2lag2,
    
    TOCC_lag1, TOCC_avg2lag1,
    #TOCC_lag1, TOCC_lag2, TOCC_avg2lag2
  )


worst.lags <-  model.dat2 %>%
  arrange(YEAR) %>%
  mutate(
    #SAM = log(SAM),
    
    COHORT_ABUND_avg2lag2 = lag(COHORT_ABUND_avg2, 2), # hard coding to lag0, no avg
    
    LG_ABUND_avg2lag2  = lag(LG_ABUND_avg2, 2),

    
    ICE_avg2lag2  = lag(ICE_avg2, 2),
    
    TOCC_avg2lag2 = lag(TOCC_avg2, 2)
    
  ) %>%
  dplyr::select(
    YEAR, PROP_INDPREF, ALL_MAT, IND_PREF,
    
    COHORT_ABUND_avg2lag2,
    
    LG_ABUND_avg2lag2,

    ICE_avg2lag2,
    TOCC_avg2lag2)


# EXPLOITATION GAM ----
# Specimen selectivity data
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  dplyr::rename(SEL = selectivity, SIZE_5MM = size) %>%
  dplyr::filter(year != "GAM predictions")

s.gam <- mgcv::gam(SEL ~ s(SIZE_5MM, k = 5),
                   method = "REML",
                   data   = sel,
                   family = Gamma(link = "log"))

# Survey specimen data
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")

spec.dat.sel <- spec.dat
spec.dat.sel$specimen <- spec.dat.sel$specimen %>%
  dplyr::mutate(
    BIN_5MM = cut_width(SIZE_1MM,
                        width  = 5,
                        center = 2.5,
                        closed = "left",
                        dig.lab = 4),
    BIN2 = BIN_5MM
  ) %>%
  tidyr::separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  dplyr::mutate(
    LOWER    = as.numeric(sub('.', '', LOWER)),
    UPPER    = as.numeric(gsub('.$', '', UPPER)),
    SIZE_5MM = (UPPER + LOWER)/2
  ) %>%
  dplyr::mutate(
    SEL             = predict(s.gam, newdata = ., type = "response"),
    SAMPLING_FACTOR = SAMPLING_FACTOR / SEL
  )

# Survey biomass, already selectivity-adjusted (>= 101 mm males)
bioabund.indpref <- crabpack::calc_bioabund(
  crab_data = spec.dat.sel,
  species   = "SNOW",
  size_min  = 101,
  size_max  = NULL,
  sex       = "male"#,
  #shell_condition = "new_hardshell"
) %>%
  mutate(
    ABUNDANCE = ABUNDANCE / 1e6,
    BIOMASS   = BIOMASS_MT / 1000        # kt
  ) %>%
  dplyr::select(YEAR, ABUNDANCE, BIOMASS)

# Directed fishery retained + discard biomass (kt)
df.dat <- read.csv("./Maturity research/Data/opilio_directedfishery_catch.csv") %>%
  mutate(directedfish_biomass = Retained_kt + Discarded_males_kt) %>%
  dplyr::select(Year, directedfish_biomass) %>%
  rename(YEAR = Year, DF_BIOMASS = directedfish_biomass)

# Natural mortality and months between survey and fishery
M <- 0.27
months_between <- 7 # Mid-survey (July) to peak fishing (January) = 6

df.exp <- df.dat %>%
  right_join(bioabund.indpref, by = "YEAR") %>%
  mutate(
    DF_BIOMASS = case_when(
      YEAR %in% c(2022:2023) ~ 0,
      TRUE ~ DF_BIOMASS
    ),
    frac_year  = months_between / 12,
    BIOMASS_fishery = BIOMASS * exp(-M * frac_year),
    EXP_RATE   = DF_BIOMASS / BIOMASS_fishery #(survey biomass available to the fishery)
  ) %>%
  dplyr::select(YEAR, DF_BIOMASS, EXP_RATE) %>%
  na.omit() %>%
  full_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1)))%>%
  filter(YEAR >=1989) %>%
  arrange(YEAR) %>%
  mutate(EXP_RATE_avg2 = zoo::rollmean(EXP_RATE,  k = 2, fill = NA, align = "right"),
         EXP_RATE_avg2lag2 = lag(EXP_RATE_avg2, 2),
         EXP_RATE_lag2 = lag(EXP_RATE, 2)) 

df.dat <- full_join(model.dat3, df.exp)

df.dat.worstlags <- full_join(worst.lags, df.exp)


# fit model 
mod <- gam(
  cbind(IND_PREF, ALL_MAT - IND_PREF) ~ 
    s(COHORT_ABUND, k = 4)+
    s(LG_ABUND_avg2, EXP_RATE_avg2, k = 10),
  data   = df.dat,
  family = quasibinomial(link = "logit"),
  method = "REML"
)

mod2 <- gam(
  cbind(IND_PREF, ALL_MAT - IND_PREF) ~ 
    s(COHORT_ABUND, k = 4)+
    s(LG_ABUND_avg2, k = 4),
  data   = df.dat,
  family = quasibinomial(link = "logit"),
  method = "REML"
)

saveRDS(mod, "./Maturity research/Models/SNOW_malepmat101_exploitation_gam_k5.rda")

worst.mod <- gam(
  cbind(IND_PREF, ALL_MAT - IND_PREF) ~ 
    s(COHORT_ABUND_avg2lag2, k = 4)+
    s(LG_ABUND_avg2lag2, EXP_RATE_avg2lag2, k = 10),
  data   = df.dat.worstlags,
  family = quasibinomial(link = "logit"),
  method = "REML"
)

saveRDS(worst.mod, "./Maturity research/Models/SNOW_malepmat101_exploitation_worstmod_k5.rda")

gam.check(mod)
summary(mod)
acf(na.omit(mod$residuals))


# smooth on a grid (as you already had)
sm.dat <- smooth_estimates(mod) %>%
  tidyr::pivot_longer(cols = 6:ncol(.), names_to = "resp", values_to = "value") %>%
  filter(.smooth == "s(COHORT_ABUND)", !is.na(value))


dat_pr <- gratia::add_partial_residuals(
  data  = mod$model,
  model = mod
)

smooth_cols <- grep("^s\\(", names(dat_pr), value = TRUE)[1]

dat_pr_long <- lapply(smooth_cols, function(sc) {
  covar <- str_match(sc, "s\\(([^)]+)\\)")[, 2]
  tibble(
    .smooth  = sc,
    x        = dat_pr[[covar]],
    .partial = dat_pr[[sc]]
  )
}) %>%
  bind_rows()

ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(sm.dat, mapping = aes(ymin = .estimate + 2 * .se, ymax = .estimate - 2 * .se), fill = "cadetblue", alpha = 0.25)+
  geom_line(color = "cadetblue", linewidth = 1.25) +
  theme_bw()+
  ylab("Partial effect")+
  xlab("Male cohort size")+
  geom_rug(data = dat_pr_long,
           aes(x = x),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.4) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        aspect.ratio = 3/4) -> cbund

p.dat <- na.omit(df.dat)
p.dat2 <- expand.grid(COHORT_ABUND = mean(p.dat$COHORT_ABUND),
                      LG_ABUND_avg2 = seq(min(p.dat$LG_ABUND_avg2), max(p.dat$LG_ABUND_avg2), length.out = 50),
                      EXP_RATE_avg2 = seq(min(p.dat$EXP_RATE_avg2), max(p.dat$EXP_RATE_avg2), length.out = 50))


preds <- cbind(p.dat2, pp = predict(mod, p.dat2, type = "response"))

# Optional: collapse to unique combinations and count how often they occur
obs_xy <- p.dat |>
  dplyr::count(LG_ABUND_avg2, EXP_RATE_avg2, name = "n")

ggplot(preds, aes(x = LG_ABUND_avg2, y = EXP_RATE_avg2)) +
  geom_tile(aes(fill = pp)) +
  geom_contour(aes(z = pp), color = "darkgrey") +
  metR::geom_text_contour(aes(z = pp), stroke = 0.2) +
  # add observed covariate combinations
  geom_point(
    data  = obs_xy,        # or p.dat if you don’t want to aggregate
    aes(x = LG_ABUND_avg2,
        y = EXP_RATE_avg2),         # drop `size = n` if you want all points the same
    inherit.aes = FALSE,
    size = 2,
    color = "hotpink",
    alpha = 0.8
  ) +
  scale_size_continuous(name = "Number of\nobservations") +  # optional
  scale_fill_viridis_c(option = "mako", name = "Proportion mature\n≥101mm") +
  theme_bw() +
  xlab("Large male abundance (t, t-1)") +
  ylab("Exploitation rate (t, t-1)") +
  theme(
    axis.text   = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title= element_text(size = 11),
    axis.title  = element_text(size = 14),
    legend.direction = "horizontal",
    legend.position  = "bottom",
    aspect.ratio = 3/4
  ) -> int

cbund/int 

ggsave("./Maturity research/Figures/expl_pmat101_k5.png", width = 6, height = 8)



