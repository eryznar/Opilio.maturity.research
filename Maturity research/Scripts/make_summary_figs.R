# SAM
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
  lab  = c("p'=0.02*", "p'=0.88")
)

ggplot(SAM.dat, aes(YEAR, SAM)) +
  geom_line() +
  geom_point() +
  ylab("Size-at-50% maturity (mm)")+
  xlab("Year")+
  geom_smooth(method = "lm") +
  facet_wrap(~ factor(sex, levels = c("Male", "Female")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14))

ggsave("./Maturity research/Figures/SNOW_SAM_timeseries.png", width = 8, height = 7)


# TARGET_SIZES ----
prop.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_INDPREF) %>%
                   rename(value = PMAT_INDPREF) %>%
                   mutate(type = "Male proportion mature ≥101mm"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, PMAT_5565) %>%
                   rename(value = PMAT_5565) %>%
                   mutate(type = "Female proportion mature 55-65mm")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, 
                           type = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm")))


m.dat <- prop.dat %>% filter(type == "Male proportion mature ≥101mm")
summary(lme(value ~ YEAR, data = na.omit(m.dat), random = ~ 1 | YEAR, correlation = corAR1()))

f.dat <- prop.dat %>% filter(type =="Female proportion mature 55-65mm")
summary(lme(value ~ YEAR, data = na.omit(f.dat), random = ~ 1 | YEAR, correlation = corAR1()))

ann_df <- data.frame(
  type   = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm"),
  x    = c(2015, 2015),  
  y    = c(0.3, 0.5),
  lab  = c("p'=0.04*", "p'=0.02*")
)

ggplot(prop.dat, aes(YEAR, value)) +
  geom_line() +
  geom_point() +
  ylab("Proportion mature")+
  xlab("Year")+
  geom_smooth(method = "lm") +
  facet_wrap(~ factor(type, levels = c("Male proportion mature ≥101mm", "Female proportion mature 55-65mm")),
             scales = "free_y", nrow = 2) +
  geom_text(data = ann_df,
            aes(x = x, y = y, label = lab),
            size = 6) +
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) 

ggsave("./Maturity research/Figures/SNOW_PMAT_timeseries.png", width = 8, height = 7)        


# MALE SAM DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm.rda")

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
  geom_histogram(breaks = seq(-15, 15, by = 5), colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 
#+ plot_annotation(title = "Male SAM diagnostics")

ggsave("./Maturity research/Figures/SNOW_maleSAM_diagnostics.png", width = 8, height = 7)

# MALE PMAT 101 DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_malepmat101_gam.rda")

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
  geom_histogram(breaks = seq(-10, 10, by = 4), colour = "black", fill = "grey80") +
  labs(title = "Histogram of residuals",
       x = "Residuals", y = "Count") +
  theme_bw()


(p_acf | p_qq) /
  (p_resid_lin | p_hist) 
#+ plot_annotation(title = "Male SAM diagnostics")

ggsave("./Maturity research/Figures/SNOW_malepmat101_diagnostics.png", width = 8, height = 7)


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

