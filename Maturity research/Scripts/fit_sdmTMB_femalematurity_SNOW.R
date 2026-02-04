## ------------------------------------------------------------
# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")


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

