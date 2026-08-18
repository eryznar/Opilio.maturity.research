## ------------------------------------------------------------
# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")

spec.dat <- crabpack::get_specimen_data(species = "SNOW", region = "EBS", channel = "API")


mod.dat <- spec.dat$specimen %>% 
              filter(SEX == 2, SHELL_CONDITION == 2) %>%
              mutate(MATURE  = case_when(CLUTCH_SIZE >0 ~ 1,
                                         TRUE ~ 0),
                     MATURE = case_when(SIZE >=72 ~ 1,
                                        SIZE <=33 ~ 0,
                                        TRUE ~ MATURE),
                     YEAR_F = as.factor(YEAR),
                     YEAR_SCALED = scale(YEAR)) %>%
          mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
                 BIN2 = BIN_5MM) %>%
          separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
          filter(SIZE <=72 & SIZE >=33) %>%
          mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                 UPPER = as.numeric(gsub('.$', '', UPPER)),
                 SIZE_5MM = (UPPER + LOWER)/2) %>%
            st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
            st_transform(., crs = "+proj=utm +zone=2") %>%
            cbind(st_coordinates(.)) %>%
            as.data.frame(.) %>%
            mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
                   LONGITUDE = X/1000)

ggplot() +
  geom_bar(
    data = mod.dat,
    mapping = aes(x = SIZE, fill = as.factor(MATURE)),
    #stat = "count",
    alpha = 0.5
  ) +
  theme_bw()+
  scale_x_continuous(breaks = seq(min(mod.dat$SIZE_5MM), max(mod.dat$SIZE_5MM), by = 10))

# Make mesh
# Set params
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")

# First evaluate between spatially varying parameterizations

# Cross validation for models with different k-smooths ----
kk<- 4:15

# Define folds
set.seed(1)
n_folds <- 5
mod.dat$fold <- sample.int(n_folds, nrow(mod.dat), replace = TRUE)


mods <- list()
results <- data.frame()
for(ii in 1:length(kk)){
  
  print(paste0("Fitting k=", kk[ii]))
  mod <- sdmTMB(MATURE ~ s(SIZE, k = kk[ii]) + YEAR_F, 
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mat.msh,
                family = binomial(),
                time = "YEAR",
                extra_time = c(2020),
                anisotropy = TRUE,
                data = mod.dat)
  
  mods[[ii]] <- mod
  
  names(mods)[ii] <- paste0("k=", kk[ii])
  
  print(paste0("CVing k=", kk[ii]))
  cv <- sdmTMB_cv(
    MATURE ~ s(SIZE, k = kk[ii]) + YEAR_F,
    data = mod.dat,
    mesh = mat.msh,
    family = binomial(),
    spatial = "on",
    spatiotemporal = "iid",
    time = "YEAR",
    extra_time = c(2020),
    anisotropy = TRUE,
    k_folds = n_folds,
    fold_ids = mod.dat$fold
  )
  
  ll <- sum(cv$fold_loglik)
  
  results <- rbind(results, data.frame(k = kk[ii], AICc(mod), logLik = ll))
  
}

results %>%
  arrange(., -logLik) -> CV.results

write.csv(CV.results, "./Maturity research/Output/femalematurity_sdmTMB_ksmooth_CV.csv")


# Cross validation for models with different mesh knots ----
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")

mod.200 <-  sdmTMB(MATURE ~ s(SIZE, k = 13) + YEAR_F, 
                          spatial = "on",
                          spatiotemporal = "iid",
                          mesh = mat.msh,
                          family = binomial(),
                          time = "YEAR",
                          extra_time = c(2020),
                          anisotropy = TRUE,
                          data = mod.dat)

cv.200 <- sdmTMB_cv(
  MATURE ~ s(SIZE, k =13) + YEAR_F,
  data = mod.dat,
  mesh = mat.msh,
  family = binomial(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "YEAR",
  extra_time = c(2020),
  anisotropy = TRUE,
  k_folds = n_folds,
  fold_ids = mod.dat$fold
)

mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

mod.300 <-  sdmTMB(MATURE ~ s(SIZE, k = 13) + YEAR_F, 
                   spatial = "on",
                   spatiotemporal = "iid",
                   mesh = mat.msh,
                   family = binomial(),
                   time = "YEAR",
                   extra_time = c(2020),
                   anisotropy = TRUE,
                   data = mod.dat)

cv.300 <- sdmTMB_cv(
  MATURE ~ s(SIZE, k =13) + YEAR_F,
  data = mod.dat,
  mesh = mat.msh,
  family = binomial(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "YEAR",
  extra_time = c(2020),
  anisotropy = TRUE,
  k_folds = n_folds,
  fold_ids = mod.dat$fold
)

# Cross validation for models with different parameterizations ----
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

mod.300spvar <-  sdmTMB(MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
                   spatial = "on",
                   spatiotemporal = "iid",
                   mesh = mat.msh,
                   family = binomial(),
                   time = "YEAR",
                   spatial_varying = ~ 0 + SIZE,
                   extra_time = c(2020),
                   anisotropy = TRUE,
                   data = mod.dat)

cv.300.spvar <- sdmTMB_cv(
  MATURE ~ s(SIZE, k =13) + YEAR_F,
  data = mod.dat,
  mesh = mat.msh,
  family = binomial(),
  spatial = "on",
  spatiotemporal = "iid",
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE_5MM,
  extra_time = c(2020),
  anisotropy = TRUE,
  k_folds = n_folds,
  fold_ids = mod.dat$fold
)

mod.300tvar <-  sdmTMB(MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
                        spatial = "on",
                        spatiotemporal = "iid",
                        mesh = mat.msh,
                        family = binomial(),
                        time = "YEAR",
                        time_varying = ~ 0 + SIZE,
                        extra_time = c(2020),
                        anisotropy = TRUE,
                        data = mod.dat)

# Best model ----
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

mod.300spvar <-  sdmTMB(MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
                        spatial = "on",
                        spatiotemporal = "iid",
                        mesh = mat.msh,
                        family = binomial(),
                        time = "YEAR",
                        spatial_varying = ~ 0 + SIZE,
                        extra_time = c(2020),
                        anisotropy = TRUE,
                        data = mod.dat)

saveRDS(mod.300spvar, "./Maturity research/Models/snowfemale_sdmTMB_spVAR_k300.rda")

# DIAGNOSTICS ----
# function
plot.resids <- function(model, model_name){
  resids <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  dat <- cbind(model$data, DHARMa_resid = resids$scaledResiduals)
  
  rr_yr  <- dat %>%
    group_by(YEAR) %>%
    arrange(DHARMa_resid, .by_group = TRUE) %>%
    mutate(
      n = n(),
      expected = ppoints(n),         # uniform quantiles
      observed = sort(DHARMa_resid)  # sort residuals for QQ
    ) %>%
    ungroup() %>%
    mutate(model = model_name)
  
  #  QQ plot with ggplot2
  ggplot()+
    theme_bw()+
    geom_point(rr_yr, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("observed")+
    xlab("expected")+
    facet_wrap(~YEAR)+
    scale_x_continuous(breaks = c(0, 0.5, 1))+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    ggtitle(model_name) -> by_yr
  
  rr_size <- dat %>%
    group_by(SIZE_5MM) %>%
    arrange(DHARMa_resid, .by_group = TRUE) %>%
    mutate(
      n = n(),
      expected = ppoints(n),         # uniform quantiles
      observed = sort(DHARMa_resid)  # sort residuals for QQ
    ) %>%
    ungroup() %>%
    mutate(model = model_name)
  
  #  QQ plot with ggplot2
  ggplot()+
    theme_bw()+
    geom_point(rr_size, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("observed")+
    xlab("expected")+
    facet_wrap(~SIZE_5MM)+
    scale_x_continuous(breaks = c(0, 0.5, 1))+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    ggtitle(model_name) -> by_size
  
  dat2 <- dat %>%
    group_by(STATION_ID) %>%
    mutate(LONGITUDE = mean(LONGITUDE), LATITUDE = mean(LATITUDE)) %>%
    ungroup()
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~YEAR)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 10)) +
    ggtitle(model_name) -> by_yr_sp
  
  ggsave(paste0("./Maturity research/Figures/snowfemale_", model_name, "spatialDHARMa_byYEAR.png"), width = 10, height = 9)
  
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~SIZE_5MM)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 10)) +
    ggtitle(model_name)-> by_size_sp
  
  ggsave(paste0("./Maturity research/Figures/snowfemale_", model_name, "spatialDHARMa_bySIZE.png"), width = 10, height = 9)
  
  return(list(by_yr = by_yr, by_size = by_size, by_yr_sp = by_yr_sp, by_size_sp = by_size_sp,
              rr_yr = rr_yr, rr_size = rr_size))
}


# Run function
model <- readRDS("./Maturity research/Models/snowfemale_sdmTMB_spVAR_k300.rda")

model_name <- "snowfem_spVAR_k300"

plot.resids(mod.300spvar, model_name) -> out
