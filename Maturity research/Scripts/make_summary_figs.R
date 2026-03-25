# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# SAM ----
SAM.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                  dplyr::select(YEAR, SAM) %>%
                   mutate(sex = "Male"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   filter(YEAR >=1989) %>%
                   mutate(sex = "Female")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, sex = c("Male", "Female")))

m.dat <- SAM.dat %>% filter(sex == "Male")
summary(lme(SAM ~ YEAR, data = na.omit(m.dat), random = ~ 1 | YEAR, correlation = corAR1()))

f.dat <- SAM.dat %>% filter(sex == "Female")
summary(lme(SAM ~ YEAR, data = na.omit(f.dat), random = ~ 1 | YEAR, correlation = corAR1()))

ann_df <- data.frame(
  sex  = c("Male", "Female"),
  x    = c(2000, 2000),   # choose positions you like
  y    = c(75, 40),
  lab  = c("p=0.02*", "p=0.88")
)

ggplot(SAM.dat, aes(YEAR, SAM)) +
  geom_line() +
  geom_point() +
  ylab("Size-at-50% maturity (mm)")+
  xlab("Year")+
  geom_smooth(method = "lm", fill = "cadetblue", color = "cadetblue", alpha = 0.2) +
  facet_wrap(~ factor(sex, levels = c("Male", "Female")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  #ggtitle("Size-at-50% maturity")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) -> SAM.plot

#ggsave("./Maturity research/Figures/SNOW_SAM_timeseries.png", width = 8, height = 7)


# TARGET_SIZES ----
prop.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_INDPREF) %>%
                   rename(value = PMAT_INDPREF) %>%
                   mutate(type = "Males ≥101mm"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_5565) %>%
                   rename(value = PMAT_5565) %>%
                   mutate(type = "Females 55-65mm")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, 
                           type = c("Males ≥101mm", "Females 55-65mm")))


m.dat <- prop.dat %>% filter(type == "Males ≥101mm")
summary(lme(value ~ YEAR, data = na.omit(m.dat), random = ~ 1 | YEAR, correlation = corAR1()))

f.dat <- prop.dat %>% filter(type =="Females 55-65mm")
summary(lme(value ~ YEAR, data = na.omit(f.dat), random = ~ 1 | YEAR, correlation = corAR1()))

ann_df <- data.frame(
  type   = c("Males ≥101mm", "Females 55-65mm"),
  x    = c(2015, 2015),  
  y    = c(0.3, 0.5),
  lab  = c("p=0.04*", "p=0.02*")
)

ggplot(prop.dat, aes(YEAR, value)) +
  geom_line() +
  geom_point() +
  ylab("Proportion mature")+
  xlab("Year")+
  geom_smooth(method = "lm", fill = "cadetblue", color = "cadetblue", alpha = 0.2) +
  facet_wrap(~ factor(type, levels = c("Males ≥101mm", "Females 55-65mm")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  #ggtitle("Proportion mature")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) -> pmat.plot

#ggsave("./Maturity research/Figures/SNOW_PMAT_timeseries.png", width = 8, height = 7)        

SAM.plot +plot_annotation(tag_levels = "A") + pmat.plot +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 17, face = "bold"),
    plot.tag.position = c(0.06, 0.98))

ggsave("./Maturity research/Figures/SNOW_PMAT&SAM_timeseries.png", width = 8, height = 6)        


# TOCC TIME SERIES FOR SI ----
tocc <- rbind(read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
            dplyr::select(YEAR, FEM_TOCC) %>%
            rename(TOCC = FEM_TOCC) %>%
            mutate(sex = "Female"),
            read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
              dplyr::select(YEAR, TOCC) %>%
              mutate(sex = "Male"))


ggplot(tocc, aes(YEAR, TOCC, color = sex))+
  geom_point(size = 2)+
  geom_line(linedwidth = 1)+
  theme_bw()+
  scale_color_manual(values = c("#8566B1", "#0095AF"), name = "")+
  ylab("Temperature occupied (°C)")+
  xlab("Year")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Maturity research/Figures/tocc_ts.png", width= 6, height = 4)



# MALE SAM DIAGNOSTICS AND EFFECT PLOTS----
mod <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm.rda")

mod2  <- mod$gam          # the GAM component
null_dev <- mod2$null.deviance
res_dev  <- mod2$deviance

dev_expl <- (null_dev - res_dev) / null_dev * 100
dev_expl

concurvity(mod$gam)

resid_df <- data.frame(
  resid   = residuals(mod$gam),                 # GAM residuals
  fitted  = fitted(mod$gam),                    # GAM fitted values
  linpred = predict(mod$gam, type = "link")     # linear predictor
)

# compute acf without plotting
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

# CI like base acf (white-noise assumption)
n   <- acf_obj$n.used
ci  <- qnorm((1 + 0.95)/2) / sqrt(n)

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = ci,  colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +  # vertical bars
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

# QQ plot
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q") +
  theme_bw()

# Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

# Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(breaks = seq(-15, 15, by = 5), colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 
#+ plot_annotation(title = "Male SAM diagnostics")

ggsave("./Maturity research/Figures/SNOW_maleSAM_diagnostics.png", width = 8, height = 7)

# Effect plots
# Add partial residuals to the original data for all smooths
dat_pr <- gratia::add_partial_residuals(
  data  = mod$gam$model,
  model = mod$gam
) 

# identify the smooth columns automatically
smooth_cols <- grep("^s\\(", names(dat_pr), value = TRUE)

dat_pr_long <- lapply(smooth_cols, function(sc) {
  # grab the covariate name from the smooth label
  # e.g. "s(COHORT_ABUND)" -> "COHORT_ABUND"
  covar <- str_match(sc, "s\\(([^)]+)\\)")[, 2]
  
  tibble(
    .smooth  = sc,
    x        = dat_pr[[covar]],
    .partial = dat_pr[[sc]]
  )
}) %>%
  bind_rows()


ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate + 1.96 * .se,
                  ymax = .estimate - 1.96 * .se),
              fill = "cadetblue", alpha = 0.25) +
  geom_line(color = "cadetblue", linewidth = 1.25) +
  # rug for all smooths
  geom_rug(data = dat_pr_long,
           aes(x = x),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.4) +
  # partial residuals for all smooths
  geom_point(data = dat_pr_long,
             aes(x = x, y = .partial),
             inherit.aes = FALSE,
             size = 1, alpha = 0.6) +
  facet_wrap(
    ~ .smooth,
    scales = "free_x",
    nrow = 2,
    labeller = as_labeller(c(
      "s(COHORT_ABUND)"  = "Male cohort abundance",
      "s(LG_ABUND_avg2)" = "Large male abundance (2-year avg)"
    ))
  ) +
  theme_bw() +
  ylab("Partial effect") +
  xlab("Value") +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )

ggsave("./Maturity research/Figures/SNOW_male_SAM_effectplots.png", width =8, height = 6)


# MALE PMAT 101 DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_malepmat101_gam.rda")

concurvity(mod)
mod$scale # should be > 1 to justify quasibinomial


# compute acf without plotting
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

# CI like base acf (white-noise assumption)
n   <- acf_obj$n.used
ci  <- qnorm((1 + 0.95)/2) / sqrt(n)

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = ci,  colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +  # vertical bars
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

# QQ plot
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q") +
  theme_bw()

# Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

# Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(bins = 7, colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 
#+ plot_annotation(title = "Male SAM diagnostics")

ggsave("./Maturity research/Figures/SNOW_malepmat101_diagnostics.png", width = 8, height = 7)

# Effect plots
# 1) Add partial residuals on link scale
dat_pr <- gratia::add_partial_residuals(
  data  = mod$model,
  model = mod,
  type  = "link",
  partial_residuals = TRUE
)

# 2) Identify smooth columns automatically
smooth_cols <- grep("^s\\(", names(dat_pr), value = TRUE)

# 3) Pivot to long: one row per obs per smooth
pr_long <- dat_pr |>
  pivot_longer(
    cols = all_of(smooth_cols),
    names_to  = ".smooth",
    values_to = ".partial"
  ) |>
  mutate(
    covariate = case_when(
      .smooth == "s(COHORT_ABUND)"  ~ COHORT_ABUND,
      .smooth == "s(LG_ABUND_avg2)" ~ LG_ABUND_avg2
      # add more cases here if you add more smooth terms
    )
  )

ggplot(sm.dat, aes(x = value, y = .estimate)) +
  geom_ribbon(aes(ymin = .estimate - 1.96 * .se,
                  ymax = .estimate + 1.96 * .se),
              fill = "cadetblue", alpha = 0.25) +
  geom_line(color = "cadetblue", linewidth = 1.25) +
  geom_rug(data = pr_long,
           aes(x = covariate),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.4) +
  geom_point(data = pr_long,
             aes(x = covariate, y = .partial),
             inherit.aes = FALSE,
             size = 1, alpha = 0.6) +
  facet_wrap(
    ~ .smooth,
    scales = "free_x",
    nrow = 2,
    labeller = as_labeller(c(
      "s(COHORT_ABUND)"  = "Male cohort abundance",
      "s(LG_ABUND_avg2)" = "Large male abundance (2-year avg)"
    ))
  ) +
  geom_smooth(data = pr_long,
              aes(x = covariate, y = .partial),
              inherit.aes = FALSE,
              method = "loess", se = FALSE, color = "grey50")+
  theme_bw() +
  xlab("Value") +
  ylab("Partial effect (link)") +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )




ggsave("./Maturity research/Figures/SNOW_male_indpref_effectplots.png", width =8, height = 6)



# FEMALE SAM DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_femaleSAM_gamm.rda")

concurvity(mod$gam)

resid_df <- data.frame(
  resid   = residuals(mod$gam),                 # GAM residuals
  fitted  = fitted(mod$gam),                    # GAM fitted values
  linpred = predict(mod$gam, type = "link")     # linear predictor
)

# compute acf without plotting
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

# CI like base acf (white-noise assumption)
n   <- acf_obj$n.used
ci  <- qnorm((1 + 0.95)/2) / sqrt(n)

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = ci,  colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +  # vertical bars
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

# QQ plot
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q") +
  theme_bw()

# Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

# Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(bins = 7, colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 

ggsave("./Maturity research/Figures/SNOW_femaleSAM_diagnostics.png", width = 8, height = 7)


# FEMALE PMAT 5565 DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_femalepmat5565_gam.rda")

concurvity(mod)
mod$scale # should be > 1 to justify quasibinomial


resid_df <- data.frame(
  resid   = residuals(mod),                 # GAM residuals
  fitted  = fitted(mod),                    # GAM fitted values
  linpred = predict(mod, type = "link")     # linear predictor
)

# compute acf without plotting
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

# CI like base acf (white-noise assumption)
n   <- acf_obj$n.used
ci  <- qnorm((1 + 0.95)/2) / sqrt(n)

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept = ci,  colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +  # vertical bars
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

# QQ plot
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Normal Q-Q") +
  theme_bw()

# Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

# Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(bins = 7, colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 

ggsave("./Maturity research/Figures/SNOW_femalepmat5565_diagnostics.png", width = 8, height = 7)


# EBS map
region_layers <- akgfmaps::get_base_layers("sebs")

## Trim survey grid to survey area
region_layers$survey.grid %>%
  st_transform(crs = st_crs(region_layers$survey.area)) %>%
  st_intersection(region_layers$survey.area) -> survey.grid

## Specify plot boundary, transform to map crs
data.frame(x = c(-178.5, -150), 
           y = c(54.5, 67)) %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(., crs = region_layers$crs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame() -> plot.boundary


## Plot 
ggplot() +
  geom_sf(data = region_layers$bathymetry, color=alpha("grey70")) +
  # geom_sf(data = region_layers$survey.grid, fill=NA, color=alpha("grey70"), linewidth = 1)+
  geom_sf(data = region_layers$survey.area, fill = alpha("darkgoldenrod", alpha=0.3), size = 0) +
  geom_sf(data = region_layers$akland, fill = "grey80", size=0.1) +
  geom_sf(data = region_layers$survey.area, fill = NA) +
  scale_x_continuous(breaks = c(-175, -165, -155)) +
  scale_y_continuous(breaks = c(56, 60, 64, 68)) +
  coord_sf(xlim = plot.boundary$X,
           ylim = plot.boundary$Y)+
  theme_bw()+
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA, color = "grey70"),
        legend.key.size = unit(0.65,'cm'),
        legend.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank()) -> study_site


# Exploitation rate ----

# Specimen selectivity data
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  dplyr::rename(SEL = selectivity, SIZE_5MM = size) %>%
  dplyr::filter(year != "GAM predictions")

s.gam <- mgcv::gam(SEL ~ s(SIZE_5MM),
                   data   = sel,
                   family = Gamma(link = "log"))

# Survey specimen data
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")

spec.dat.sel <- spec.dat
spec.dat.sel$specimen <- spec.dat.sel$specimen %>%
  dplyr::mutate(
    BIN_5MM = ggplot2::cut_width(SIZE_1MM,
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
  sex       = "male"
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
  filter(YEAR >=1989)
 
summary(lme(EXP_RATE ~ YEAR, data = na.omit(df.exp), random = ~ 1 | YEAR, correlation = corAR1()))

ggplot(df.exp, aes(YEAR, EXP_RATE))+
  geom_point()+
  geom_line()+
  theme_bw()+
  #geom_smooth(method = "lm")+
  ylab("Exploitation rate")+
  xlab("Year")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) -> exp.rate


# Ice ----
# Load Jan-April ice data
ice <- read.csv(paste0("./Maturity research/Output/ebs_ice_means_1980-", current.year, ".csv")) %>%
  dplyr::select(year, se, value) %>%
  rename(YEAR = year, ICE = value) %>%
  filter(YEAR > 1988)

ggplot(ice, aes(YEAR, ICE))+
  geom_point()+
  geom_line()+
  #geom_errorbar(aes(ymin = ICE - se, ymax = ICE + se))+
  theme_bw()+
  xlab("Year")+
  ylab("Ice area fraction")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) -> ice.plot

# Male abundance plots ----
ind.pref <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                     size_min = 101, size_max = NULL,  sex = "male",
                                     shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND) %>%
  mutate(cat = "Industry-preferred males (≥101mm)")


lg.male <- crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                   size_min = 95, size_max = NULL,  sex = "male") %>%
                                   #shell_condition = c("new_hardshell", "oldshell", "very_oldshell")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)) %>% # convert to kt
  dplyr::select(YEAR, ABUND) %>%
  mutate(cat = "Large males (≥95mm)")

cohort <- crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                   size_min = 40, size_max = 60,  sex = "male", 
                                  shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)) %>% # convert to kt
  dplyr::select(YEAR, ABUND) %>%
  mutate(cat = "Pre-mature males (40-60mm)")

rbind(cohort, lg.male, ind.pref) %>%
  group_by(cat) %>%
  mutate(ABUND_SCALED = scale(ABUND)) %>%
  ungroup() %>%
  filter(YEAR > 1988) %>%
  full_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), cat = unique(.$cat))) -> plot.dat


ggplot(plot.dat, aes(YEAR, ABUND_SCALED, color = cat))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  theme_bw()+
  scale_color_manual(values = c("#26185F","#0095AF", "#9ADCBB"), name = "")+
  xlab("Year")+
  ylab("Scaled abundance (millions)")+
  theme(
    legend.position = c(0.3, 0.99),   # x, y in npc (0–1)
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text= element_text(size =10)) -> male.plot 

# Female abundance plots ----
cohort <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 35, size_max = 45,  sex = "female", 
                                    shell_condition = c("new_hardshell")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND) %>%
  mutate(cat = "Pre-mature females (35-45mm)")


imm_3545 <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                   size_min = 35, size_max = 45,  sex = "female", 
                                   crab_category = c("immature_female")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND)

imm <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                     size_min = NULL, size_max = NULL,  sex = "female", 
                                     crab_category = c("immature_female")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND)

all <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                size_min = NULL, size_max = NULL,  sex = "female") %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND)

pimm<- imm$ABUND/all$ABUND
pimm3545<- imm_3545$ABUND/all$ABUND

cor(cohort$ABUND, pimm) # cohort strength vs prop immature
cor(cohort$ABUND, pimm3545) # cohort strength vs prop immature 3545
cor(cohort$ABUND, all$ABUND) # cohort strength vs all females



# mature female abundance
matfem <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                                size_min = NULL, size_max = NULL,  sex = "female", 
                                                crab_category = c("mature_female")) %>%
  group_by(YEAR) %>%
  reframe(ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, ABUND) %>%
  mutate(cat = "Mature females")


rbind(cohort, matfem) %>%
  group_by(cat) %>%
  mutate(ABUND_SCALED = scale(ABUND)) %>%
  ungroup() %>%
  filter(YEAR > 1988) %>%
  full_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), cat = unique(.$cat))) -> plot.dat


ggplot(plot.dat, aes(YEAR, ABUND_SCALED, color = cat))+
  geom_point(size = 2)+
  geom_line(linewidth = 1)+
  theme_bw()+
  scale_color_manual(values = c("#8566B1", "#ABC4DE"), name = "")+
  xlab("Year")+
  ylab("Scaled abundance (millions)")+
  theme(
    legend.position = c(0.2, 0.99),   # x, y in npc (0–1)
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = NA),
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.text= element_text(size =10)) -> female.plot 

# Big fig 1
#study site, ice time series, and exploitation rate time series. Then one panel plotting male abundance 
# (pre-maturity / 40-60mm, >95, and >101) - all of these scaled so highest abundance year = 1. 
# Then one panel for female abundance (35-45 mm and mature) also scaled to 1.  
(male.plot/female.plot + plot_layout(axis_titles = "collect", axes = "collect_x"))
ptop <- (study_site + ice.plot + exp.rate)+  
        theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "lines"),
              plot.tag.position = c(0.25, 0.98),   # x, y in [0,1]
              plot.tag = element_text(face = "bold", size = 11))
                       

pbottom <- male.plot/female.plot + plot_layout(axis_titles = "collect", axes = "collect_x", heights = c(1.4, 1.4)) +
  theme(plot.margin = unit(c(0, 0.2, 0.2, 0.2), "lines"),
        plot.tag.position = c(0.25, 0.98),   # x, y in [0,1]
        plot.tag = element_text(face = "bold", size = 11))

ptop/pbottom + plot_layout(heights = c(1, 2)) 


ggsave("./Maturity research/Figures/Fig1.png", height= 6, width = 8)
