
# MALE SAM DIAGNOSTICS ----
mod <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm.rda")

## 1. Concurvity of smooths (GAM part)
concurvity(mod$gam, full = TRUE)   # or 'para = TRUE' etc. as needed
# or use gam.check(mod$gam) / qq.gam(mod$gam, rep = 100) for extra QQ bands

## 2. Residual diagnostics for the GAMM (lme part)

# use normalized residuals from the lme component
resid_df <- data.frame(
  resid   = residuals(mod$lme, type = "normalized"),
  fitted  = fitted(mod$lme),
  linpred = predict(mod$lme, type = "link")
)

## ACF of normalized residuals
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

n  <- acf_obj$n.used
ci <- qnorm((1 + 0.95)/2) / sqrt(n)   # white-noise CI, matches base acf()

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept =  ci, colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

## QQ plot of normalized residuals
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q") +
  theme_bw()

## Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

## Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(bins = 7, colour = "black", fill = "grey80") +
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

## 1. Concurvity of smooths (GAM part)
concurvity(mod$gam, full = TRUE)   # or 'para = TRUE' etc. as needed
# or use gam.check(mod$gam) / qq.gam(mod$gam, rep = 100) for extra QQ bands

## 2. Residual diagnostics for the GAMM (lme part)

# use normalized residuals from the lme component
resid_df <- data.frame(
  resid   = residuals(mod$lme, type = "normalized"),
  fitted  = fitted(mod$lme),
  linpred = predict(mod$lme, type = "link")
)

## ACF of normalized residuals
acf_obj <- acf(resid_df$resid, plot = FALSE)
acf_df  <- with(acf_obj, data.frame(lag = lag, acf = acf))

n  <- acf_obj$n.used
ci <- qnorm((1 + 0.95)/2) / sqrt(n)   # white-noise CI, matches base acf()

p_acf <- ggplot(acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, colour = "black") +
  geom_hline(yintercept =  ci, colour = "blue", linetype = "dashed") +
  geom_hline(yintercept = -ci, colour = "blue", linetype = "dashed") +
  geom_segment(aes(xend = lag, y = 0, yend = acf)) +
  labs(title = "ACF",
       x = "Lag", y = "ACF") +
  theme_bw()

## QQ plot of normalized residuals
p_qq <- ggplot(resid_df, aes(sample = resid)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  labs(title = "Q-Q") +
  theme_bw()

## Residuals vs linear predictor
p_resid_lin <- ggplot(resid_df, aes(x = linpred, y = resid)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "red") +
  labs(title = "Residuals vs linear predictor",
       x = "Linear predictor", y = "Residuals") +
  theme_bw()

## Histogram of residuals
p_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(bins = 6, colour = "black", fill = "grey80") +
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

