## ------------------------------------------------------------
# PURPOSE: to generate a simulation to assess drivers of declining proportion mature at 101mm

# Author: Emily Ryznar


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD MODEL DATA/MODEL ------------------------------------------------------
m.dat <- read.csv("./Maturity research/Output/indpref_modeldat.csv") %>%
          dplyr::select(!X) 


# DATA INPUTS ----
# Specimen data and selectivity
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  rename(SEL = selectivity, SIZE_5MM = size) %>%
  filter(year != "GAM predictions")

s.gam <- gam(SEL ~ s(SIZE_5MM), data = sel, family = Gamma(link = "log"))

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

# Statistical model (pmat_101 GAM)
mod <- readRDS("./Maturity research/Models/SNOW_malepmat101_gam.rda")

# Ogives
ogive_raw <- read.csv("./Maturity research/Data/SNOW_maleogives_withselectivity.csv")

# Historical 40-60mm survey abundance by 1mm size bins
cohort <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                    size_min = 40, size_max = 60,  sex = "male", 
                                    shell_condition = c("new_hardshell"), bin_1mm = TRUE) %>%
  group_by(YEAR, SIZE_1MM) %>%
  reframe(COHORT_ABUND = sum(ABUNDANCE)/1e6) %>% # convert to kt
  dplyr::select(YEAR, SIZE_1MM, COHORT_ABUND)  %>%
  filter(YEAR >= 1989)

# Exploitation rate
bioabund.indpref <-  crabpack::calc_bioabund(crab_data = spec.dat.sel, species = "SNOW", 
                                             size_min = 101, size_max = NULL,  sex = "male") %>%
  mutate(ABUNDANCE = ABUNDANCE/1e6,
         BIOMASS = BIOMASS_MT/1000) %>% # convert to kt
  dplyr::select(YEAR, ABUNDANCE, BIOMASS)

df.dat <- read.csv("./Maturity research/Data/opilio_directedfishery_catch.csv") %>%
  mutate(directedfish_biomass = Retained_kt+ Discarded_males_kt) %>% #
  dplyr::select(Year, directedfish_biomass) %>%
  rename(YEAR = Year, DF_BIOMASS = directedfish_biomass) %>%
  right_join(., bioabund.indpref) %>% # to calculate exploitation rate
  mutate(DF_BIOMASS = case_when((YEAR %in% c(2020, 2022:2023)) ~ NA,
                                TRUE ~ DF_BIOMASS),
         EXP_RATE = DF_BIOMASS/BIOMASS) %>%
  dplyr::select(YEAR, DF_BIOMASS, EXP_RATE) %>%
  na.omit()

## PREPARE OGIVE MATRIX ----
ogive <- ogive_raw %>%
  dplyr::select(YEAR, SIZE_5MM, PROP_MATURE) %>%
  rename(p_term = PROP_MATURE)

years_vec <- sort(unique(ogive$YEAR))
sizes_5   <- sort(unique(ogive$SIZE_5MM))

ogive_mat <- matrix(NA_real_,
                    nrow = length(years_vec),
                    ncol = length(sizes_5),
                    dimnames = list(years_vec, sizes_5))

for (i in seq_len(nrow(ogive))) {
  y  <- ogive$YEAR[i]
  sz <- ogive$SIZE_5MM[i]
  ogive_mat[as.character(y), as.character(sz)] <- ogive$p_term[i]
}

## HELPER FUNCTIONS ----
# use the ogive range, but keep everything inside it
sizes_5   <- sort(unique(ogive$SIZE_5MM))          # e.g., 7.5 ... 172.5
min_size  <- ceiling(min(sizes_5))                 # 8
max_size  <- floor(max(sizes_5))                   # 172
size_bins <- min_size:max_size                     # 8:172
n_size    <- length(size_bins)


get_p_term_year <- function(year_val, size_bins, ogive_mat) {
  ogive_years <- as.numeric(rownames(ogive_mat))
  
  # if exact year not present, return NA vector
  if (!year_val %in% ogive_years) {
    return(rep(NA_real_, length(size_bins)))
  }
  
  row_idx <- which(ogive_years == year_val)
  
  valid_5mm <- as.numeric(colnames(ogive_mat))
  col_idx <- sapply(size_bins, function(sz) which.min(abs(valid_5mm - sz)))
  
  p_vec <- ogive_mat[row_idx, col_idx]
  p_vec[is.na(p_vec)] <- 0
  p_vec
}

## GROWTH AND MORTALITY ----
a_male <- -15 # intercept from assessment growth fit
b_male <- 0.40 # positive slope from assessment

mean_growth <- function(size_mm) {
  pmax(0, a_male + b_male * size_mm)
}


M <- 0.271 # mortality from SAFE (applying HAmel's methodology)
surv_M <- exp(-M)

## GAM-BASED P(mature at 101 mm) ----
p_pref_101 <- function(cohort_abund, lg_abund_avg2, mod) {
  newdat <- data.frame(
    COHORT_ABUND  = cohort_abund,
    LG_ABUND_avg2 = lg_abund_avg2
  )
  predict(mod, newdata = newdat, type = "response")[1]
}

## RECRUITMENT FROM COHORT DF ----
# cohort: YEAR, SIZE_1MM, COHORT_ABUND (long format)
years_cohort <- sort(unique(cohort$YEAR))

draw_recruits <- function(cohort) {
  yr <- sample(unique(cohort$YEAR), size = 1)
  sub <- cohort %>% filter(YEAR == yr)
  rec_vec <- sub$COHORT_ABUND
  names(rec_vec) <- sub$SIZE_1MM
  rec_vec
}

## ONE-YEAR PROJECTION FUNCTION ----
project_one_year <- function(N_t,
                             year_val,
                             expl_rate,
                             mod,
                             size_bins,
                             ogive_mat,
                             M,
                             cohort) {
  
  n_size <- length(size_bins)
  
  N_surv <- N_t * exp(-M)
  
  inc <- mean_growth(size_bins)
  new_size <- size_bins + inc
  new_idx <- findInterval(new_size, vec = size_bins, all.inside = TRUE)
  
  N_grown <- numeric(n_size)
  for (i in seq_along(size_bins)) {
    N_grown[new_idx[i]] <- N_grown[new_idx[i]] + N_surv[i]
  }
  
  large_idx   <- which(size_bins >= 101)
  catch_large <- N_grown[large_idx] * expl_rate
  N_post_fish <- N_grown
  N_post_fish[large_idx] <- N_grown[large_idx] - catch_large
  
  p_term <- get_p_term_year(year_val, size_bins, ogive_mat)
  
  if (all(is.na(p_term))) {
    # no ogive this year: carry N_forward but mark outputs as NA
    term_molts <- rep(0, length(size_bins))
    N_after_tm <- N_post_fish
    cohort_abund_t  <- sum(N_post_fish[size_bins >= 40 & size_bins <= 60])
    lg_abund_avg2_t <- sum(N_post_fish[size_bins >= 95])
    
    return(list(
      N_next          = N_after_tm,
      term_molts      = term_molts,
      catch_large     = catch_large,
      cohort_abund_t  = cohort_abund_t,
      lg_abund_avg2_t = lg_abund_avg2_t,
      has_ogive       = FALSE
    ))
  }
  
  idx_101 <- which(size_bins == 101)
  if (length(idx_101) == 1) {
    cohort_abund_t  <- sum(N_post_fish[size_bins >= 40 & size_bins <= 60])
    lg_abund_avg2_t <- sum(N_post_fish[size_bins >= 95])
    p_mature_101_t  <- p_pref_101(cohort_abund_t, lg_abund_avg2_t, mod)
    p_term[idx_101] <- p_mature_101_t
  }
  
  term_molts <- N_post_fish * p_term
  N_after_tm <- N_post_fish - term_molts
  
  rec <- draw_recruits(cohort)
  rec_sizes <- as.numeric(names(rec))
  idx_rec   <- match(rec_sizes, size_bins)
  keep      <- !is.na(idx_rec)
  N_after_tm[idx_rec[keep]] <- N_after_tm[idx_rec[keep]] + rec[keep]
  
  # state variables for output (after fishing, before next year)
  cohort_abund_t  <- sum(N_post_fish[size_bins >= 40 & size_bins <= 60])
  lg_abund_avg2_t <- sum(N_post_fish[size_bins >= 95])
  
  list(
    N_next          = N_after_tm,
    term_molts      = term_molts,
    catch_large     = catch_large,
    cohort_abund_t  = cohort_abund_t,
    lg_abund_avg2_t = lg_abund_avg2_t
  )
}

# SIMULATION FUNCTION ----
simulate_scenario <- function(expl_rate,
                              years_sim,
                              start_year,
                              mod,
                              size_bins,
                              ogive_mat,
                              cohort,
                              M,
                              replicate_id = 1) {
  
  n_size <- length(size_bins)
  
  N0 <- numeric(n_size)
  rec0 <- draw_recruits(cohort)
  rec_sizes0 <- as.numeric(names(rec0))
  idx0 <- match(rec_sizes0, size_bins)
  keep0 <- !is.na(idx0)
  N0[idx0[keep0]] <- rec0[keep0]
  
  N_t <- N0
  
  prop_tm_above101 <- numeric(years_sim)
  cohort_abund_ts  <- numeric(years_sim)
  lg_abund_avg2_ts <- numeric(years_sim)
  sim_years        <- numeric(years_sim)
  
  for (t in 1:years_sim) {
    sim_year <- start_year + t - 1
    
    res <- project_one_year(
      N_t         = N_t,
      year_val    = sim_year,
      expl_rate   = expl_rate,
      mod         = mod,
      size_bins   = size_bins,
      ogive_mat   = ogive_mat,
      M           = M,
      cohort      = cohort
    )
    
    tm <- res$term_molts
    if (!is.null(res$has_ogive) && !res$has_ogive) {
      prop_tm_above101[t] <- NA
    } else {
      idx_large_tm <- which(size_bins >= 101)
      total_tm <- sum(tm)
      tm_above <- sum(tm[idx_large_tm])
      prop_tm_above101[t] <- if (total_tm > 0) tm_above / total_tm else NA
    }
    idx_large_tm <- which(size_bins >= 101)
    total_tm <- sum(tm)
    tm_above <- sum(tm[idx_large_tm])
    prop_tm_above101[t] <- if (total_tm > 0) tm_above / total_tm else NA
    
    cohort_abund_ts[t]  <- res$cohort_abund_t
    lg_abund_avg2_ts[t] <- res$lg_abund_avg2_t
    sim_years[t]        <- sim_year
    
    N_t <- res$N_next
  }
  
  list(
    replicate_id      = replicate_id,
    year              = sim_years,
    prop_tm_above101  = prop_tm_above101,
    cohort_abund_101  = cohort_abund_ts,
    lg_abund_avg2     = lg_abund_avg2_ts
  )
}

## RUN EXAMPLE SCENARIOS ----
years_sim  <- length(years_vec)+5
start_year <- min(years_vec)
end_year <- max(years_vec)
expl_vec   <- seq(0, 1, by = 0.2)
n_reps <- 100  # 100 independent simulations per exploitation rate

sim_df <- do.call(rbind, lapply(seq_along(expl_vec), function(i) {
  U <- expl_vec[i]
  
  # run n_reps replicates for this exploitation rate
  reps <- lapply(1:n_reps, function(r) 
    simulate_scenario(U,
                      years_sim,
                      start_year,
                      mod,
                      size_bins,
                      ogive_mat,
                      cohort,
                      M,
                      replicate_id = r)
  )
  
  # stack replicates into a data.frame with replicate column
  do.call(rbind, lapply(reps, function(res) {
    data.frame(
      exploitation     = U,
      replicate        = res$replicate_id,
      year             = res$year,
      prop_tm_above101 = res$prop_tm_above101,
      cohort_abund_101 = res$cohort_abund_101,
      lg_abund_avg2    = res$lg_abund_avg2
    )
  }))
}))

sim_df %>%
  filter(prop_tm_above101 >0) %>%
  group_by(exploitation, year) %>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(exploitation, year) %>%
  reframe(pind = mean(prop_tm_above101),
          N = n(),
          se = sd(prop_tm_above101)/sqrt(N)) %>%
  full_join(., expand.grid(year = seq(min(.$year), max(.$year)), exploitation = expl_vec)) -> sum_df


ggplot(sum_df %>% filter(exploitation !=1), aes(year, pind, color = as.factor(exploitation)))+
  geom_line()+
  geom_point()+
  geom_errorbar(sum_df%>% filter(exploitation !=1), mapping = aes(year, ymin = pind-se, ymax = pind+se, color = as.factor(exploitation)))+
  theme_bw()



ggplot(sim_df %>% filter(prop_tm_above101>0), aes(cohort_abund_101, lg_abund_avg2, fill = prop_tm_above101))+
  geom_point(shape = 21, stroke = NA, size =2, alpha = 0.5)+
  facet_wrap(~exploitation)+
  scale_fill_viridis_c()+
  theme_bw()
