## ------------------------------------------------------------
## PURPOSE: Simulation to assess drivers of declining proportion
##          mature (terminal molt) for snow crab males >=101 mm
##          using time steps, with GAM vs ogive comparison.
## Author: Emily Ryznar
## ------------------------------------------------------------

## LOAD LIBRARIES / PARAMS ---------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

## LOAD MODEL DATA / MODELS --------------------------------------------------

m.dat <- read.csv("./Maturity research/Output/indpref_modeldat.csv") %>%
  dplyr::select(!X)

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

# GAM for proportion of newly mature crab that are >=101 mm
mod <- readRDS("./Maturity research/Models/SNOW_malepmat101_gam.rda")

# Ogives (historical, by year)
ogive_raw <- read.csv("./Maturity research/Data/SNOW_maleogives_withselectivity.csv")

# Historical 40–60 mm abundance (for recruitment sampling)
cohort <- crabpack::calc_bioabund(
  crab_data       = spec.dat.sel,
  species         = "SNOW",
  size_min        = 40,
  size_max        = 60,
  sex             = "male",
  shell_condition = c("new_hardshell"),
  bin_1mm         = TRUE
) %>%
  dplyr::group_by(YEAR, SIZE_1MM) %>%
  dplyr::reframe(COHORT_ABUND = sum(ABUNDANCE) / 1e6) %>%  # kt
  dplyr::select(YEAR, SIZE_1MM, COHORT_ABUND) %>%
  dplyr::filter(YEAR >= 1989)

## PREPARE OGIVE MATRIX ------------------------------------------------------
ogive <- ogive_raw %>%
  dplyr::select(YEAR, SIZE_5MM, PROP_MATURE) %>%
  dplyr::rename(p_term = PROP_MATURE)

## MEAN OGIVE ACROSS YEARS (DEFAULT SHAPE) --------------------------------
mean_ogive_df <- ogive %>%
  dplyr::group_by(SIZE_5MM) %>%
  dplyr::summarise(p_term_mean = mean(p_term), .groups = "drop")

get_p_term_mean <- function(size_bins, mean_ogive_df) {
  valid_5mm <- mean_ogive_df$SIZE_5MM
  p5        <- mean_ogive_df$p_term_mean
  col_idx   <- sapply(size_bins, function(sz) which.min(abs(valid_5mm - sz)))
  p_vec     <- p5[col_idx]
  p_vec[is.na(p_vec)] <- 0
  p_vec
}

sizes_5   <- sort(unique(ogive$SIZE_5MM))
min_size  <- ceiling(min(sizes_5))
max_size  <- floor(max(sizes_5))
size_bins <- min_size:max_size
n_size    <- length(size_bins)

p_term_mean_all <- get_p_term_mean(size_bins, mean_ogive_df)

years_vec <- sort(unique(ogive$YEAR))

ogive_mat <- matrix(
  NA_real_,
  nrow = length(years_vec),
  ncol = length(sizes_5),
  dimnames = list(as.character(years_vec), as.character(sizes_5))
)

for (i in seq_len(nrow(ogive))) {
  y  <- ogive$YEAR[i]
  sz <- ogive$SIZE_5MM[i]
  ogive_mat[as.character(y), as.character(sz)] <- ogive$p_term[i]
}

get_p_term_year <- function(year_label, size_bins, ogive_mat) {
  ogive_years <- rownames(ogive_mat)
  if (!as.character(year_label) %in% ogive_years) {
    return(rep(NA_real_, length(size_bins)))
  }
  row_idx   <- which(ogive_years == as.character(year_label))
  valid_5mm <- as.numeric(colnames(ogive_mat))
  col_idx   <- sapply(size_bins, function(sz) which.min(abs(valid_5mm - sz)))
  p_vec     <- ogive_mat[row_idx, col_idx]
  p_vec[is.na(p_vec)] <- 0
  p_vec
}

## GROWTH AND NATURAL MORTALITY ---------------------------------------------
a_male <- -15
b_male <- 0.40

mean_growth <- function(size_mm) {
  pmax(0, a_male + b_male * size_mm)
}

M      <- 0.271
surv_M <- exp(-M)

## GAM-BASED PROPORTION (NEWLY MATURE >=101) --------------------------------
p_pref_new_ge101 <- function(cohort_abund, lg_abund_avg2, mod) {
  newdat <- data.frame(
    COHORT_ABUND  = cohort_abund,
    LG_ABUND_avg2 = lg_abund_avg2
  )
  predict(mod, newdata = newdat, type = "response")[1]
}

## RECRUITMENT DRAW ----------------------------------------------------------
draw_recruits <- function(cohort) {
  yr  <- sample(unique(cohort$YEAR), size = 1)
  sub <- cohort %>% dplyr::filter(YEAR == yr)
  rec_vec <- sub$COHORT_ABUND
  names(rec_vec) <- sub$SIZE_1MM
  rec_vec
}

## IMPLIED PROPORTION (NEWLY MATURE >=101) FOR AN OGIVE YEAR ----------------
## Ogive-implied newly mature >=101 proportion = (new mature >=101)/(all new mature)
implied_new_mature_ge101 <- function(ogive_year,
                                     N_vec,
                                     size_bins,
                                     ogive_mat) {
  p_term_y <- get_p_term_year(ogive_year, size_bins, ogive_mat)
  N_mature_by_size <- N_vec * p_term_y
  total_mature     <- sum(N_mature_by_size)
  if (total_mature <= 0) return(NA_real_)
  idx_ge101      <- which(size_bins >= 101)
  mature_ge101   <- sum(N_mature_by_size[idx_ge101])
  mature_ge101 / total_mature
}

## CHOOSE OGIVE YEAR MATCHING GAM TARGET (NEWLY MATURE SCALE) ---------------
choose_ogive_year_from_N <- function(target_prop_new_ge101,
                                     N_vec,
                                     size_bins,
                                     ogive_mat) {
  ogive_years <- as.numeric(rownames(ogive_mat))
  if (length(ogive_years) == 0L) return(NA_real_)
  props <- sapply(ogive_years, function(yr) {
    implied_new_mature_ge101(yr, N_vec, size_bins, ogive_mat)
  })
  if (all(is.na(props))) return(NA_real_)
  diffs <- abs(props - target_prop_new_ge101)
  ogive_years[which.min(diffs)]
}

## ONE-TIME-STEP PROJECTION  -------------------------------
## Order: mortality -> exploitation (>=101) -> growth -> (GAM/ogive) -> terminal molt -> recruitment
project_one_step <- function(N_t,        # immature / not-yet-terminal
                             M_t,        # already mature
                             expl_rate,
                             mod,
                             size_bins,
                             ogive_mat,
                             mean_ogive_df,
                             M,
                             cohort) {
  
  n_size <- length(size_bins)
  
  # 1. natural mortality
  N_surv <- N_t * exp(-M)
  M_surv <- M_t * exp(-M)
  
  # 2. exploitation on survivors (>=101 mm), before growth
  idx_fished  <- which(size_bins >= 101)
  catch_imm   <- N_surv[idx_fished] * expl_rate
  catch_mat   <- M_surv[idx_fished] * expl_rate
  
  N_post_fish0 <- N_surv
  M_post_fish0 <- M_surv
  N_post_fish0[idx_fished] <- N_surv[idx_fished] - catch_imm
  M_post_fish0[idx_fished] <- M_surv[idx_fished] - catch_mat
  
  # 3. growth (same kernel)
  inc      <- mean_growth(size_bins)
  new_size <- size_bins + inc
  new_idx  <- findInterval(new_size, vec = size_bins, all.inside = TRUE)
  
  N_grown <- numeric(n_size)
  M_grown <- numeric(n_size)
  for (i in seq_along(size_bins)) {
    N_grown[new_idx[i]] <- N_grown[new_idx[i]] + N_post_fish0[i]
    M_grown[new_idx[i]] <- M_grown[new_idx[i]] + M_post_fish0[i]
  }
  
  # After growth, these are the post-fishing abundances by size
  N_post_fish <- N_grown
  M_post_fish <- M_grown
  
  # 4. state variables for GAM (total abundance)
  N_total_post <- N_post_fish + M_post_fish
  cohort_abund_t  <- sum(N_total_post[size_bins >= 40 & size_bins <= 60])
  lg_abund_avg2_t <- sum(N_total_post[size_bins >= 95])
  
  # 5. GAM target = newly mature >=101 proportion
  target_prop_new_ge101 <- p_pref_new_ge101(cohort_abund_t, lg_abund_avg2_t, mod)
  
  # 6. choose ogive year (newly mature scale)
  ogive_year_star <- choose_ogive_year_from_N(
    target_prop_new_ge101,
    N_vec     = N_total_post,
    size_bins = size_bins,
    ogive_mat = ogive_mat
  )
  
  # 7. get p_term and ogive-implied newly mature proportion
  if (is.na(ogive_year_star)) {
    p_term <- get_p_term_mean(size_bins, mean_ogive_df)
    N_mature_new_by_size <- N_total_post * p_term
    total_mature_new     <- sum(N_mature_new_by_size)
    if (total_mature_new <= 0) {
      prop_new_ge101_ogive <- NA_real_
    } else {
      idx_ge101            <- which(size_bins >= 101)
      mature_ge101_new     <- sum(N_mature_new_by_size[idx_ge101])
      prop_new_ge101_ogive <- mature_ge101_new / total_mature_new
    }
  } else {
    p_term <- get_p_term_year(ogive_year_star, size_bins, ogive_mat)
    prop_new_ge101_ogive <- implied_new_mature_ge101(
      ogive_year_star, N_total_post, size_bins, ogive_mat
    )
  }
  
  # 8. apply terminal molt: move from N_post_fish to M_post_fish
  term_molts <- N_post_fish * p_term          # newly mature this year, by size
  N_after_tm <- N_post_fish - term_molts
  M_after_tm <- M_post_fish + term_molts      # includes newly mature
  
  # 8b. newly mature >=101 proportion from the simulated term_molts
  total_new_mature_step <- sum(term_molts)
  if (total_new_mature_step > 0) {
    idx_ge101_step        <- which(size_bins >= 101)
    new_mature_ge101_step <- sum(term_molts[idx_ge101_step])
    prop_new_ge101_sim    <- new_mature_ge101_step / total_new_mature_step
  } else {
    prop_new_ge101_sim <- NA_real_
  }
  
  # 9. add recruits into immature state
  rec       <- draw_recruits(cohort)
  rec_sizes <- as.numeric(names(rec))
  idx_rec   <- match(rec_sizes, size_bins)
  keep      <- !is.na(idx_rec)
  N_after_tm[idx_rec[keep]] <- N_after_tm[idx_rec[keep]] + rec[keep]
  
  # 10. stock-level proportion mature >=101 (M state)
  total_mature_stock <- sum(M_after_tm)
  if (total_mature_stock > 0) {
    idx_ge101_stock    <- which(size_bins >= 101)
    mature_ge101_stock <- sum(M_after_tm[idx_ge101_stock])
    prop_ge101_stock   <- mature_ge101_stock / total_mature_stock
  } else {
    prop_ge101_stock <- NA_real_
  }
  
  list(
    N_next                      = N_after_tm,
    M_next                      = M_after_tm,
    term_molts                  = term_molts,
    catch_large_imm             = catch_imm,
    catch_large_mat             = catch_mat,
    cohort_abund_t              = cohort_abund_t,
    lg_abund_avg2_t             = lg_abund_avg2_t,
    prop_new_mature_ge101_gam   = target_prop_new_ge101,   # GAM (newly mature)
    prop_new_mature_ge101_ogive = prop_new_ge101_ogive,    # ogive (newly mature)
    prop_new_mature_ge101_sim   = prop_new_ge101_sim,      # simulated newly mature
    prop_mature_ge101_stock     = prop_ge101_stock,        # stock-level
    ogive_year_star             = ogive_year_star
  )
}


## SIMULATION FUNCTION (TIME-STEP BASED) ------------------------------------
## Takes N0, M0 so we can share N0 across exploitation within replicate
simulate_scenario <- function(expl_rate,
                              n_steps,
                              mod,
                              size_bins,
                              ogive_mat,
                              mean_ogive_df,
                              cohort,
                              M,
                              replicate_id = 1,
                              burn_in = 0,
                              N0,
                              M0 = NULL) {
  
  n_size <- length(size_bins)
  if (is.null(M0)) M0 <- numeric(n_size)
  
  message(paste("Replicate", replicate_id, "- U =", expl_rate))
  
  N_t <- N0
  M_t <- M0
  
  time_step                  <- 1:n_steps
  prop_tm_above101           <- numeric(n_steps)  # stock-scale
  cohort_abund_ts            <- numeric(n_steps)
  lg_abund_avg2_ts           <- numeric(n_steps)
  prop_new_ge101_gam_ts      <- numeric(n_steps)
  prop_new_ge101_ogive_ts    <- numeric(n_steps)
  prop_new_ge101_sim_ts      <- numeric(n_steps)
  prop_ge101_stock_ts        <- numeric(n_steps)
  ogive_year_star_ts         <- numeric(n_steps)
  
  for (t in 1:n_steps) {
    
    res <- project_one_step(
      N_t          = N_t,
      M_t          = M_t,
      expl_rate    = expl_rate,
      mod          = mod,
      size_bins    = size_bins,
      ogive_mat    = ogive_mat,
      mean_ogive_df = mean_ogive_df,
      M            = M,
      cohort       = cohort
    )
    
    prop_tm_above101[t]        <- res$prop_mature_ge101_stock
    cohort_abund_ts[t]         <- res$cohort_abund_t
    lg_abund_avg2_ts[t]        <- res$lg_abund_avg2_t
    prop_new_ge101_gam_ts[t]   <- res$prop_new_mature_ge101_gam
    prop_new_ge101_ogive_ts[t] <- res$prop_new_mature_ge101_ogive
    prop_new_ge101_sim_ts[t]   <- res$prop_new_mature_ge101_sim
    prop_ge101_stock_ts[t]     <- res$prop_mature_ge101_stock
    ogive_year_star_ts[t]      <- res$ogive_year_star
    
    N_t <- res$N_next
    M_t <- res$M_next
  }
  
  out <- data.frame(
    replicate_id               = replicate_id,
    time_step                  = time_step,
    exploitation               = expl_rate,
    # stock-scale output
    prop_tm_above101           = prop_tm_above101,          # mature stock ≥101 / all mature
    cohort_abund_101           = cohort_abund_ts,
    lg_abund_avg2              = lg_abund_avg2_ts,
    # newly-mature scale outputs
    prop_new_mature_ge101_gam   = prop_new_ge101_gam_ts,
    prop_new_mature_ge101_ogive = prop_new_ge101_ogive_ts,
    prop_new_mature_ge101_sim   = prop_new_ge101_sim_ts,
    # stock-scale again for clarity
    prop_mature_ge101_stock     = prop_ge101_stock_ts,
    ogive_year_star             = ogive_year_star_ts
  )
  
  # drop burn-in time steps (e.g., burn_in = 15)
  if (burn_in > 0) {
    out <- out[out$time_step > burn_in, ]
  }
  
  out
}

## RUN SCENARIOS -------------------------------------------------------------
n_steps <- 50
expl_vec <- seq(0, 0.8, by = 0.2)
n_reps   <- 1000
burn_in  <- 5

sim_list <- vector("list", n_reps)

for (r in 1:n_reps) {
  
  n_size <- length(size_bins)
  
  # initial immature from recruitment (replicate-specific)
  N0_r <- numeric(n_size)
  rec0 <- draw_recruits(cohort)
  rec_sizes0 <- as.numeric(names(rec0))
  idx0 <- match(rec_sizes0, size_bins)
  keep0 <- !is.na(idx0)
  N0_r[idx0[keep0]] <- rec0[keep0]
  
  # initial mature stock for this replicate
  M0_r <- numeric(n_size)
  
  # run all exploitation scenarios with this N0_r, M0_r
  sim_list[[r]] <- do.call(rbind, lapply(expl_vec, function(U)
    simulate_scenario(
      expl_rate     = U,
      n_steps       = n_steps,
      mod           = mod,
      size_bins     = size_bins,
      ogive_mat     = ogive_mat,
      mean_ogive_df = mean_ogive_df,
      cohort        = cohort,
      M             = M,
      replicate_id  = r,
      burn_in       = burn_in,
      N0            = N0_r,
      M0            = M0_r
    )
  ))
}

sim_df <- do.call(rbind, sim_list)
write.csv(sim_df, "./Maturity research/Output/indpref_sim_output.csv")

## SUMMARY FOR TIME-SERIES PLOT ---------------------------------------------
# Cumulative proportion mature ≥101
sum_df <- sim_df %>%
  dplyr::filter(!is.na(prop_tm_above101)) %>%
  dplyr::group_by(exploitation, time_step) %>%
  dplyr::summarise(
    pind_med = stats::median(prop_tm_above101),
    pind_lo  = stats::quantile(prop_tm_above101, 0.1, na.rm = TRUE),
    pind_hi  = stats::quantile(prop_tm_above101, 0.90, na.rm = TRUE),
    N        = dplyr::n(),
    .groups  = "drop"
  )

ggplot(sum_df %>% dplyr::filter(exploitation != 1),
       aes(time_step, pind_med, color = as.factor(exploitation))) +
  geom_ribbon(aes(ymin = pind_lo, ymax = pind_hi,
                  fill = as.factor(exploitation)),
              alpha = 0.2, color = NA) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ggtitle("Cumulative proportion mature ≥101mm")+
  labs(x = "Time step (year)",
       y = "Proportion mature stock ≥101",
       color = "Exploitation",
       fill  = "Exploitation")

# Newly matured
new_sum_df <- sim_df %>%
  dplyr::filter(!is.na(prop_new_mature_ge101_sim)) %>%
  dplyr::group_by(exploitation, time_step) %>%
  dplyr::reframe(
    pind_med = stats::median(prop_new_mature_ge101_sim),
    pind_lo  = stats::quantile(prop_new_mature_ge101_sim, 0.025, na.rm = TRUE),
    pind_hi  = stats::quantile(prop_new_mature_ge101_sim, 0.975, na.rm = TRUE),
    N        = dplyr::n(),
    se   = stats::sd(prop_new_mature_ge101_sim) / sqrt(N),
    pind_mean = mean(prop_new_mature_ge101_sim)
  )

ggplot(new_sum_df %>% dplyr::filter(exploitation != 1),
       aes(time_step, pind_med, color = as.factor(exploitation))) +
  geom_ribbon(aes(ymin = pind_lo, ymax = pind_hi,
                  fill = as.factor(exploitation)),
              alpha = 0.2, color = NA) +
  geom_line() +
  geom_point() +
  #facet_wrap(~exploitation)+
  theme_bw() +
  ggtitle("Newly mature entering ≥101mm")+
  labs(x = "Time step (year)",
       y = "Proportion newly mature ≥101",
       color = "Exploitation",
       fill  = "Exploitation")

ggplot(new_sum_df %>% dplyr::filter(exploitation != 1, time_step > 15),
       aes(time_step, pind_mean, color = as.factor(exploitation))) +
  geom_errorbar(aes(ymin = pind_mean - se, ymax = pind_mean + se,
                color = as.factor(exploitation))) +
  geom_line() +
  geom_point() +
  #facet_wrap(~exploitation)+
  theme_bw() +
  ggtitle("Newly mature entering ≥101mm")+
  labs(x = "Time step (year)",
       y = "Proportion newly mature ≥101",
       color = "Exploitation",
       fill  = "Exploitation")

## STATE SPACE PLOTS ---------------------------------------------------------
ggplot(sim_df %>% dplyr::filter(!is.na(prop_tm_above101), exploitation !=1),
       aes(cohort_abund_101, lg_abund_avg2, fill = prop_tm_above101)) +
  geom_point(shape = 21, stroke = NA, size = 2, alpha = 0.5) +
  facet_wrap(~exploitation) +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Cumulative proportion mature ≥101mm")+
  labs(x = "40–60 mm cohort abundance",
       y = "≥95 mm abundance",
       fill = "Prop mature ≥101mm")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggplot(sim_df %>%
         dplyr::filter(!is.na(prop_new_mature_ge101_sim),
                       exploitation != 1),
       aes(cohort_abund_101, lg_abund_avg2,
           fill = prop_new_mature_ge101_sim)) +
  geom_point(shape = 21, stroke = NA, size = 2, alpha = 0.5) +
  facet_wrap(~exploitation) +
  scale_fill_viridis_c() +
  ggtitle("Newly mature entering ≥101mm")+
  theme_bw() +
  labs(x = "40–60 mm cohort abundance",
       y = "≥95 mm abundance",
       fill = "Prop newly mature ≥101")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

## GAM vs OGIVE (NEWLY MATURE SCALE) ----------------------------------------
gam_ogive_df <- sim_df %>%
  dplyr::filter(!is.na(prop_new_mature_ge101_gam),
                !is.na(prop_new_mature_ge101_ogive),
                prop_new_mature_ge101_ogive > 0) %>%
  dplyr::group_by(exploitation, time_step) %>%
  dplyr::summarise(
    prop_ge101_gam_mean   = mean(prop_new_mature_ge101_gam,   na.rm = TRUE),
    prop_ge101_ogive_mean = mean(prop_new_mature_ge101_ogive, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(gam_ogive_df %>% dplyr::filter(exploitation < 1),
       aes(prop_ge101_gam_mean, prop_ge101_ogive_mean)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE, method = "lm") +
  theme_bw() +
  facet_wrap(~exploitation, scales = "free")+
  labs(x = "GAM predicted proportion newly mature ≥101",
       y = "Ogive-implied proportion newly mature ≥101")
