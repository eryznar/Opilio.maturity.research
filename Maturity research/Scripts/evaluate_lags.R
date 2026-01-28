model.dat <- read.csv("./Maturity research/Data/model_dat.csv") %>%
  dplyr::select(!c(X, EXP_RATE)) 

#1) build covariate variants
max_lag <- 6


dat <- model.dat %>%
  arrange(YEAR) %>%
  mutate(
    ICE_avg2 = rollmean(ICE, k = 2, fill = NA, align = "right"),
    ICE_avg3 = rollmean(ICE, k = 3, fill = NA, align = "right"),

    INST1_ABUND_avg2 = rollmean(INST1_ABUND, k = 2, fill = NA, align = "right"),
    INST1_ABUND_avg3 = rollmean(INST1_ABUND, k = 3, fill = NA, align = "right"),

    
    LG_ABUND_avg2  = rollmean(LG_ABUND, k = 2, fill = NA, align = "right"),
    LG_ABUND_avg3  = rollmean(LG_ABUND, k = 3, fill = NA, align = "right"),

    TOCC_avg2 = rollmean(TOCC, k = 2, fill = NA, align = "right"),
    TOCC_avg3 = rollmean(TOCC, k = 3, fill = NA, align = "right")
  )

#2) Extract ccf metrics for each series
vars <- names(dat)[-c(1,2)]
cc_df <- data.frame()

for (ii in seq_along(vars)) {
  pp <- dat[[vars[ii]]]
  
  ok <- complete.cases(dat$SAM, pp)
  if (sum(ok) < 2) next   # skip if not enough data
  
  SAM_ok <- dat$SAM[ok]
  pp_ok  <- pp[ok]
  
  cc <- ccf(SAM_ok, pp_ok, lag.max = max_lag, plot = FALSE)
  
  df <- data.frame(
    var = vars[ii],
    lag = as.numeric(cc$lag),
    cor = as.numeric(cc$acf)
  )
  
  cc_df <- rbind(cc_df, df)
}

#3) specify running means
long_df <- cc_df %>%
              mutate(smooth = case_when(grepl("avg2", var) ~ "2-year",
                                        grepl("avg3", var) ~ "3-year",
                                        TRUE ~ "none"),
                     var = case_when(grepl("avg2", var) ~ gsub("_avg2", "", var),
                                     grepl("avg3", var) ~ gsub("_avg3", "", var),
                                     TRUE ~ var))


#4) plot
ggplot(long_df, aes(lag, cor, fill = factor(smooth, levels = c("none", "2-year", "3-year"))))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("cadetblue", "salmon", "darkgoldenrod"), name = "smooth")+
  facet_wrap(~var)+
  theme_bw()+
  scale_x_continuous(breaks = seq(-max_lag, max_lag, by = 1), labels = seq(-max_lag, max_lag, 1))+
  theme(
  panel.grid.minor.x =element_blank())

#4) select best avg-lag combo for each covariate
long_df %>%
  filter(lag >= 0) %>%
  mutate(cor_abs = abs(cor)) %>%
  group_by(var) %>%
  slice_max(order_by = cor_abs, n = 3, with_ties = FALSE) %>%
  ungroup()

#5) add in important lags
dat2 <- dat %>% 
  mutate(LG_ABUND_lag3 = lag(LG_ABUND, 3),
         LG_ABUND_avg3lag5 = lag(LG_ABUND_avg3, 5))

#5) run models
# Fit models
response <- "SAM"

# Generate parameter combos to fit over
lg.pars <-  c(NA, "LG_ABUND_lag3", "LG_ABUND", "LG_ABUND_avg3lag5")
sm.pars <- c(NA, "INST1_ABUND_avg2", "INST1_ABUND", "INST1_ABUND_avg3")
tocc.pars <-  c(NA, "TOCC_avg3", "TOCC", "TOCC_avg3")
ice.pars <- c(NA, "ICE", "ICE_avg2", "ICE_avg3")

combos <- expand_grid(
  #df = df.pars,
  ice = ice.pars,
  lg = lg.pars,
  sm = sm.pars,
  tocc = tocc.pars) %>%
  filter(!(is.na(lg) & is.na(sm) & is.na(ice) & is.na(tocc)))

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
    fit <- gamm(fml, data = dat2, family = gaussian(), 
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

mod <- gamm(SAM ~ s(INST1_ABUND_avg3, k =8)+
              s(LG_ABUND_avg3lag5, k = 4)+
              s(TOCC_avg3, k = 4)+
              s(ICE_avg3, k = 4),
            correlation = corAR1(),
            data = dat2,
            family = gaussian())

mod <- gam(SAM ~ s(INST1_ABUND_avg2, k = 4)+
              s(LG_ABUND, k = 4)+
              s(TOCC_avg3, k = 4),
            data = dat2,
            family = gaussian())

diagnose.gamm(mod)
diagnose(mod)
