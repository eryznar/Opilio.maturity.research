# PURPOSE: To evaluate the proportion of sh2 Chionoecetes mature males that are targeted by the fishery and to evaluate
# depletion of large males by the fishery for Chionoecetes and RKC

# Author: Emily Ryznar

# LOAD LIBRARIES ----
library(tidyverse)
library(crabpack)
library(patchwork)
library(nlme)
library(lme4)
library(mgcv)

# GET SPECIMEN DATA ----
# Set channel
channel <- "KOD"
years <- 1989:2025

# Pull specimen data
species <- "SNOW"
snow_spec <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = years, # set in libs/params script
                                             channel = channel)


# SNOW CRAB Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Calculate biomass/abundance 
snow.ind.pref <- calc_bioabund(crab_data = snow_spec,
                          species = "SNOW",
                          years = years, # set in libs/params script
                          sex = "male",
                          shell_condition = "new_hardshell",
                          size_min = 101) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = snow_spec,
                         species = "SNOW",
                         years = years, # set in libs/params script
                         sex = "male",
                         shell_condition = "new_hardshell")


# Add in mature male abundance from the tech memo 2025
tm_matabund <- data.frame(tm_matmaleabund = c(458.3, 678.2, 380.5, 418.8, 666.6, 487.8, 597.8, 592.6, 
                                              286.1, 225.8, 82.1, 98.1, 209.9, 116.1, 119, 189.2, 386, 215.9, 339.8,
                                              NA, 150.9, 400.9, 267.8, NA, 255.7, NA, 178.2, NA, 454.3,
                                              1232.4, 559.2, 115.3, 82.6, 49.9, 297.8, 574.4), 
                          ci = c(104, 192, 124.7, 96.3, 211.1, 112.4, 122.9, 116.2, 66.6, 45.3, 18.6, 26.9, 74.3,
                                 61.1, 33, 57.4, 117, 73.4, 84.2, NA, 30.8, 85.1, 76.1, NA, 56.5, NA, 54, NA, 146.1,
                                 367.8, 188.8, 41, 34.7, 14, 80.8, 154.5),
                          YEAR = c(1989:2019, 2021:2025)) %>%
  mutate(tm_matmaleabund = tm_matmaleabund*1e6)

indpref.matmaleprop.snow <- snow.ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)


# TANNER CRAB WEST Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Pull specimen data
species <- "TANNER"
tan_spec <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = years, # set in libs/params script
                                             channel = channel)

# Calculate biomass/abundance 
tanW.ind.pref <- calc_bioabund(crab_data = tan_spec,
                          species = "TANNER",
                          #years = years, # set in libs/params script
                          sex = "male",
                          district = "W166",
                          shell_condition = "new_hardshell",
                          size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = tan_spec,
                         species = "TANNER",
                         years = years, # set in libs/params script
                         sex = "male",
                         district = "W166",
                         shell_condition = "new_hardshell")

# Add in mature male abundance from the tech memo 2025
tm_matabund <- data.frame(tm_matmaleabund = c(70.3, 40.4, 16.5, 21, 8.2, 4.6, 3.1, 1.8, 4.1, 5.1, 3.9, 7.9, 8.1, 14.3, 16.8, 40,
                               39.2, 38.7, 48.1, NA, 30.9, 22, 12.8, NA, 53.6, NA, 26.7, 6.9, 13, 9, 13.3,
                               10.9, 20.2, 35.4, 87.5), 
                          ci = c(32.7, 17.5, 5.7, 6.8, 3.2, 1.4, 1.3, 0.6, 1.7, 4.2, 1.5, 3.3, 3.7, 5.5, 10, 14.2, 16.3,
                                 11.8, 27.4, NA, 16.5, 19.1, 5.3, NA, 19.2, NA, 10.1, 2.4, 3.9, 2.2, 4.6, 3.2, 5, 12.6, 36.1),
                          YEAR = c(1990:2019, 2021:2025)) %>%
  mutate(tm_matmaleabund = tm_matmaleabund*1e6)


indpref.matmaleprop.tanW <- tanW.ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)

# TANNER CRAB EAST Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Calculate biomass/abundance 
tanE.ind.pref <- calc_bioabund(crab_data = tan_spec,
                          species = "TANNER",
                          years = years, # set in libs/params script
                          sex = "male",
                          district = "E166",
                          shell_condition = "new_hardshell",
                          size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = tan_spec,
                         species = "TANNER",
                         years = years, # set in libs/params script
                         sex = "male",
                         district = "E166",
                         shell_condition = "new_hardshell")

# Add in mature male abundance from the tech memo 2025
tm_matabund <- data.frame(tm_matmaleabund = c(47.3, 36.6, 47.5, 23.4, 12.6, 1.1, 1.2, 1.5, 4.7, 5.9, 8.2,
                                              5.9, 1.9, 4.8, 6.4, 10.2, 16.5, 15.2, 32.2, NA, 5.9, NA, 31.0, NA, 
                                              39.2, NA, 8.1, 2.7, 0.7, 1.9, 7, 13.4, 8, 12.5, 9.6), 
                          ci = c(17, 33.2, 42.5, 13.6, 6.3, 0.9, 0.6, 0.7, 1.8, 4.9, 7.7, 3.2, 1.1, 2.5, 2.6, 3.4, 
                                 21.2, 13.7, 41.7, NA, 4.3, NA, 20.8, NA, 15.1, NA, 5.9, 4.1, 0.4, 1.9, 3.6, 6.1, 
                                 3.2, 7.4, 3.8),
                          YEAR = c(1990:2019, 2021:2025)) %>%
  mutate(tm_matmaleabund = tm_matmaleabund*1e6)

indpref.matmaleprop.tanE <- tanE.ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)

# JOIN PLOTS ----
rbind(indpref.matmaleprop.snow %>% mutate(species = "Snow"), 
      indpref.matmaleprop.tanE %>% mutate(species = "Tanner East"), 
      indpref.matmaleprop.tanW %>% mutate(species = "Tanner West")) -> plot.dat


plot.dat <- right_join(plot.dat, data.frame(YEAR = rep(seq(min(plot.dat$YEAR), max(plot.dat$YEAR)),3),
                                            species = rep(c("Snow", "Tanner East", "Tanner West"), each = 37)))

plot.dat %>% filter(species == "Snow") -> snow
plot.dat %>% filter(species == "Tanner East") -> tanE
plot.dat %>% filter(species == "Tanner West") -> tanW

snow.gamm <- gamm(indpref_matmaleprop ~ s(YEAR, k = 4), data = snow, correlation = corAR1(form = ~YEAR))
tanE.gamm <- gamm(indpref_matmaleprop ~ s(YEAR, k = 4), data = tanE, correlation = corAR1(form = ~YEAR))
tanW.gamm <- gamm(indpref_matmaleprop ~ s(YEAR, k = 4), data = tanW, correlation = corAR1(form = ~YEAR))

summary(snow.gamm$gam)
gam.check(snow.gamm$gam)
summary(tanE.gamm$gam)
gam.check(tanE.gamm$gam)
summary(tanW.gamm$gam)
gam.check(tanW.gamm$gam)

p2 <- predict(snow.gamm$gam, snow, se = TRUE) %>%
      cbind(snow, .) %>%
      rename(pred = fit, pred.se = se.fit) %>%
      mutate(pred.ci = 1.96*pred.se)

p3 <- predict(tanE.gamm$gam, tanE, se = TRUE) %>%
        cbind(tanE, .) %>%
        rename(pred = fit, pred.se = se.fit) %>%
        mutate(pred.ci = 1.96*pred.se)

p4 <- predict(tanW.gamm$gam, tanW, se = TRUE) %>%
  cbind(tanW, .) %>%
  rename(pred = fit, pred.se = se.fit) %>%
  mutate(pred.ci = 1.96*pred.se)

ggplot() +
  geom_line(plot.dat, mapping = aes(YEAR, indpref_matmaleprop))+
  geom_point(plot.dat, mapping = aes(YEAR, indpref_matmaleprop))+
  geom_ribbon(p2 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  geom_line(p2, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  geom_ribbon(p3 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  geom_line(p3, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  geom_ribbon(p4 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  geom_line(p4, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  ylab("Proportion industry preferred")+
  facet_wrap(~species, nrow = 3, scales = "free_y")+
  theme_bw()+
  xlab("Year")+
  theme(axis.text = element_text(size= 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))

ggsave("./Figures/chionoecetes_propindpref_matmale.png", width = 7, height = 7)


plot.dat %>% dplyr::select(species, YEAR, indpref_abund, ABUNDANCE_CI, tm_matmaleabund, ci) %>%
  rename("Industry-preferred" = indpref_abund, indpref_ci = ABUNDANCE_CI, "Mature male" = tm_matmaleabund, matmale_ci = ci) -> pp 
  
right_join(pp %>%
                dplyr::select(species, YEAR, "Industry-preferred", "Mature male") %>%
              pivot_longer(., cols = c(3, 4), names_to = "type", values_to = "abundance") %>%
                mutate(abundance = abundance/1e6),
            pp %>%
              dplyr::select(species, YEAR, indpref_ci, matmale_ci) %>%
              mutate(indpref_ci = indpref_ci/1e6) %>%
              rename("Industry-preferred" = indpref_ci, "Mature male" = matmale_ci) %>%
              pivot_longer(., cols = c(3, 4), names_to = "type", values_to = "ci")) -> pd


ggplot() +
  geom_ribbon(pd , mapping = aes(x = YEAR, ymin = abundance - ci, ymax = abundance + ci, fill = type), alpha = 0.25)+
  geom_line(pd, mapping = aes(YEAR, abundance, color = type), linewidth = 1)+
  geom_point(pd, mapping = aes(YEAR, abundance, color = type))+
  ylab("Abundance (millions)")+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue"))+
  scale_fill_manual(values = c("darkgoldenrod", "cadetblue"))+
  facet_wrap(~species, nrow = 3, scales = "free_y")+
  theme_bw()+
  xlab("Year")+
  theme(axis.text = element_text(size= 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))


ggsave("./Figures/chionoecetes_indpref_matmale_timeseries.png", width = 7, height = 7)


# SNOW ind preferred biomass vs. df catch biomass ----
snow.ind.pref <- calc_bioabund(crab_data = snow_spec,
                               species = "SNOW",
                               years = years, # set in libs/params script
                               sex = "male",
                               #shell_condition = "new_hardshell",
                               size_min = 101) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

df.catch <- data.frame(YEAR = 1982:2023, catch_kt = c(11.85, 12.16, 29.94, 44.45, 46.22, 61.4, 67.79, 73.4, 149.1, 143, 104.7, 67.94,
                                                      34.13, 39.81, 54.22, 114.4, 88.09, 15.1, 11.46, 14.8, 12.84, 10.86, 11.29, 16.77, 
                                                      16.49, 28.59, 26.56, 21.78, 24.61, 40.29, 30.05, 24.49, 30.82, 18.42, 9.67, 8.6,
                                                      12.51, 15.43, 20.41, 2.48, NA, NA))

snow.mod.dat <- right_join(snow.ind.pref, df.catch)

snow.mod.dat %>%
  dplyr::select(SPECIES, YEAR, indpref_bio, catch_kt) %>%
  mutate(survey_indpref_bio_kt = indpref_bio*0.001,
          retained_catch_kt = catch_kt,
          ratio = (retained_catch_kt/survey_indpref_bio_kt) * 100) %>%
  dplyr::select(!c(indpref_bio, catch_kt)) -> snow.mod.dat2
# 
# cody.dat <- read.csv("./Data/obs_exp_dat.csv")
# 
# ggplot() +
#   geom_line(snow.mod.dat2 %>% filter(YEAR >=1989), mapping = aes(YEAR, ratio))+
#   geom_point(snow.mod.dat2 %>% filter(YEAR >=1989), mapping = aes(YEAR, ratio))+
#   geom_line(cody.dat %>% filter(currency == ">101 mm"), mapping = aes(year, (tot_retained_wt/biomass *100)), col = "blue") +
#   geom_point(cody.dat %>% filter(currency == ">101 mm"), mapping = aes(year, (tot_retained_wt/biomass *100)), col = "blue") +
#   # geom_ribbon(p3 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
#   # geom_line(p3, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
#   ylab("Retained catch/survey industry-preferred biomass")+
#   theme_bw()+
#   xlab("Year")+
#   theme(axis.text = element_text(size= 12),
#         axis.title = element_text(size = 12),
#         strip.text = element_text(size = 12))
# 


# Tanner ind preferred biomass vs. df catch biomass ----
tanE.ind.pref <- calc_bioabund(crab_data = tan_spec,
                               species = "TANNER",
                               district = "E166",
                               years = years, # set in libs/params script
                               sex = "male",
                               #shell_condition = "new_hardshell",
                               size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

tanW.ind.pref <- calc_bioabund(crab_data = tan_spec,
                               species = "TANNER",
                               district = "W166",
                               years = years, # set in libs/params script
                               sex = "male",
                               #shell_condition = "new_hardshell",
                               size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)


df.catch <- data.frame(YEAR = c(1980:1996, 2005:2023), 
                       catchE = c(rep(NA, 17), 0, 631.2, 710, 806.9, 592.4, 0, 0, 0, 643.3, 3829.3, 5107.7, 0, 
                                  0.1, 0, 0, 0, 0, 528.4, 343.5),
                       catchW = c(rep(NA, 17), 244.5, 155.5, 151.1, 47.2, 0, 0, 0, 0, 593.6, 2368.7, 3770.3, 0, 1117.5,
                                  1103.9, 0, 655.2, 493.5, 384.9, 596.8),
                       All = c(13426, 4990, 2390, 549, 1429, NA, NA, 998, 3180, 11113, 18189, 14424, 15921,
                               7666, 3538, 1919, 821, rep(NA, length(2005:2023)))) %>%
             mutate(All = case_when((YEAR %in% 2005:2023) ~ catchE+catchW,
                                    TRUE ~ All)) %>%
            pivot_longer(., cols = c(catchE, catchW, All), names_to = "species", values_to = "df_bio") %>%
            mutate(species = case_when((species == "catchE") ~ "Tanner East",
                                       (species == "catchW") ~ "Tanner West",
                                       TRUE ~ "Tanner"))
                         
                         
rbind(tanE.ind.pref %>% 
      dplyr::select(YEAR, indpref_bio) %>% 
      mutate(species = "Tanner East") %>% 
      rename(survey_indpref_bio = indpref_bio),
      tanW.ind.pref %>% 
        dplyr::select(YEAR, indpref_bio) %>% 
        mutate(species = "Tanner West") %>%
        rename(survey_indpref_bio = indpref_bio)) -> tcatch

tcatch %>%
  group_by(YEAR) %>%
  reframe(survey_indpref_bio = sum(survey_indpref_bio)) %>%
  mutate(species = "Tanner") -> allcatch

tanner.mod.dat <- rbind(tcatch, allcatch) %>%
                  right_join(., df.catch) %>%
                  mutate(prop = (df_bio/survey_indpref_bio)*100)



## Join Chionoecetes data ----
dat <- rbind(snow.mod.dat2 %>% 
               mutate(SPECIES = "Snow") %>%
               rename(survey_indpref_bio = survey_indpref_bio_kt,
                      df_bio = retained_catch_kt,
                      prop = ratio,
                      species = SPECIES) , tanner.mod.dat) %>%
      filter(YEAR >= 1989) %>%
      right_join(data.frame(YEAR = rep(seq(min(.$YEAR), max(.$YEAR)), 2), species = rep(c("Snow", "Tanner", "Tanner East", "Tanner West"), 
                                                                                        each = 35))) %>%
      mutate(prop = case_when((is.infinite(prop) == TRUE) ~ NA,
                              TRUE ~ prop))

dat %>% filter(species == "Snow") %>% na.omit() -> snow
dat %>% filter(species == "Tanner East")  %>% na.omit() -> tanE
dat %>% filter(species == "Tanner West")  %>% na.omit() -> tanW
dat %>% filter(species == "Tanner")  %>% na.omit() -> tan

snow.gamm <- gamm(prop ~ s(YEAR), data = snow, correlation = corAR1(form = ~YEAR))
tanE.gamm <- gamm(prop ~ s(YEAR), data = tanE, correlation = corAR1(form = ~YEAR))
tanW.gamm <- gamm(prop ~ s(YEAR), data = tanW, correlation = corAR1(form = ~YEAR))
tan.gamm <- gamm(prop ~ s(YEAR, k =4), data = tan, correlation = corAR1(form = ~YEAR))

summary(snow.gamm$gam)
gam.check(snow.gamm$gam)
summary(tanE.gamm$gam)
gam.check(tanE.gamm$gam)
summary(tanW.gamm$gam)
gam.check(tanW.gamm$gam)
summary(tan.gamm$gam)
gam.check(tan.gamm$gam)

p3 <- predict(snow.gamm$gam, se = TRUE) %>%
  cbind(snow, .) %>%
  rename(pred = fit, pred.se = se.fit) %>%
  mutate(pred.ci = 1.96*pred.se)

p4 <- predict(tan.gamm$gam, se = TRUE) %>%
  cbind(tan, .) %>%
  rename(pred = fit, pred.se = se.fit) %>%
  mutate(pred.ci = 1.96*pred.se)


ggplot() +
  geom_line(dat %>% filter(is.infinite(prop) == FALSE, !species %in% c("Tanner East", "Tanner West"),
                           YEAR >=1989), mapping = aes(YEAR, prop))+
  geom_point(dat %>% filter(is.infinite(prop) == FALSE, !species %in% c("Tanner East", "Tanner West"),
                            YEAR >=1989), mapping = aes(YEAR, prop))+
  geom_ribbon(p3 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  geom_line(p3, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  # geom_ribbon(p4 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  # geom_line(p4, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  ylab("Retained catch/survey industry-preferred biomass")+
  facet_wrap(~species, nrow = 2, scales = "free_y")+
  theme_bw()+
  xlab("Year")+
  theme(axis.text = element_text(size= 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))

ggsave("./Figures/chionoecetes_exploitationrate_indprefmale.png", width = 7, height = 7)


# BBRKC ind preferred biomass vs. df catch biomass ----
# # Get specimen data
# specimen_data <- crabpack::get_specimen_data(species = "RKC",
#                                              region = "EBS",
#                                              district = "BB",
#                                              years = 1979:2024, # set in libs/params script
#                                              channel = "API")
# # Calculate biomass/abundance 
# ind.pref <- calc_bioabund(crab_data = specimen_data,
#                           species = "RKC",
#                           years = 1979:2024, # set in libs/params script
#                           sex = "male",
#                           district = "BB",
#                           #shell_condition = "new_hardshell",
#                           size_min = 135) %>%
#   rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)
# 
# df.catch <- data.frame(YEAR = c(1979:2022), catch_kt = c(107.83, 129.95, 33.37, 2.99, NA, 4.08, 4.09, 11.31, 12.29, 7.36, 10.16, 20.44,
#                                                          17.18, 8.07, 14.59, NA, NA, 8.5, 8.91, 15, 11.84, 8.24, 8.52, 9.67, 15.73, 15.45,
#                                                          18.31, 15.62, 20.37, 20.33, 15.93, 14.83, 7.83, 7.85, 8.6, 9.99, 9.97, 8.47, 6.6,
#                                                          4.31, 3.79, 2.65, NA, NA))
# 
# mod.dat <- right_join(ind.pref, df.catch)
# 
# mod <- lm(log(mod.dat$indpref_bio + 10) ~ log(mod.dat$catch_kt + 10), data = mod.dat)
# 
# ggplot(mod.dat, aes(log(catch_kt + 10), log(indpref_bio + 10)))+
#   geom_point()+
#   geom_smooth(method = "lm", se = FALSE) + 
#   theme_bw()+
#   ylab("log(survey biomass + 10)")+
#   xlab("log(commercial catch biomass + 10)")+
#   ggtitle("BBRKC ≥ 135mm survey biomass vs. commercial catch biomass") -> p.6
# 
# 
# p.4/p.5


# SIZE AT 50% MATURITY TIMESERIES ----
species <- "SNOW"
snow_mat <- crabpack::get_male_maturity(species = species,
                                         region = "EBS",
                                         #years = years, # set in libs/params script
                                         channel = channel)

snow_L50 <- snow_mat$model_parameters %>%
              dplyr::select(YEAR, B_EST, B_SE)

species <- "TANNER"
tanE_mat <- crabpack::get_male_maturity(species = species,
                                        region = "EBS",
                                        district = "E166",
                                        #years = years, # set in libs/params script
                                        channel = channel)

tanE_L50 <- tanE_mat$model_parameters %>%
  dplyr::select(YEAR, B_EST, B_SE)

tanW_mat <- crabpack::get_male_maturity(species = species,
                                        region = "EBS",
                                        district = "W166",
                                        #years = years, # set in libs/params script
                                        channel = channel)

tanW_L50 <- tanW_mat$model_parameters %>%
  dplyr::select(YEAR, B_EST, B_SE)



snow_mod <- gamm(B_EST ~ s(YEAR), data = snow_L50, correlation = corAR1(form = ~YEAR))
summary(snow_mod$gam)

tanE_mod <- gamm(B_EST ~ s(YEAR), data = tanE_L50, correlation = corAR1(form = ~YEAR))
summary(tanE_mod$gam)

tanW_mod <- gamm(B_EST ~ s(YEAR), data = tanW_L50, correlation = corAR1(form = ~YEAR))
summary(tanW_mod$gam)


p3 <- cbind(predict(snow_mod$gam, snow_L50, se = TRUE), snow_L50) %>%
      rename(pred = fit, pred.se = se.fit) %>%
      mutate(pred.ci = 1.96*pred.se, species = "Snow")

L50 <- rbind(snow_L50 %>% mutate(species = "Snow"), 
             tanE_L50 %>% mutate(species = "Tanner East"), 
             tanW_L50 %>% mutate(species = "Tanner West")) %>%
       right_join(data.frame(YEAR = rep(seq(min(.$YEAR), max(.$YEAR)),3),
                             species = rep(c("Snow", "Tanner East", "Tanner West"), each = 37)))

ggplot() +
  geom_line(L50, mapping = aes(YEAR, B_EST))+
  geom_point(L50, mapping = aes(YEAR, B_EST))+
  geom_line(p3, mapping = aes(YEAR, pred), color = "cadetblue", linewidth = 1)+
  geom_errorbar(L50, mapping = aes(YEAR, ymax = B_EST + B_SE, ymin = B_EST - B_SE))+
  geom_ribbon(p3 , mapping = aes(x = YEAR, ymin = pred - pred.ci, ymax = pred + pred.ci), fill = "cadetblue", alpha = 0.25)+
  ylab("Size at 50% maturity (mm)")+
  facet_wrap(~species, nrow = 3, scales = "free_y")+
  theme_bw()+
  xlab("Year")+
  theme(axis.text = element_text(size= 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))

ggsave("./Figures/chionoecetes_l50timeseries.png", width = 7, height = 7)

## EXPLORATORY ANALYSES ----
pp <- right_join(L50, plot.dat)

p1 <- pp %>% filter(species == "Snow")
p2 <- pp %>% filter(species == "Tanner East")
p3 <- pp %>% filter(species == "Tanner West")


m1 <- gamm(B_EST ~ s(indpref_matmaleprop), na.omit(p1), correlation = corAR1(form = ~YEAR))
m2 <- gls(B_EST ~ indpref_matmaleprop, na.omit(p2), correlation = corAR1(form = ~YEAR))
m3 <- gls(B_EST ~ indpref_matmaleprop, na.omit(p3), correlation = corAR1(form = ~YEAR))

pr1 <- cbind(pred = predict(m1, na.omit(p1)), na.omit(p1)) %>% mutate(species = "Snow")
pr2 <- cbind(pred = predict(m2, na.omit(p2)), na.omit(p2)) %>% mutate(species = "Tanner East")



ggplot()+
  geom_point(pp, mapping = aes(indpref_matmaleprop, B_EST))+
  geom_line(pr1, mapping = aes(indpref_matmaleprop, pred), color = "cadetblue", linewidth = 1)+
  geom_line(pr2, mapping = aes(indpref_matmaleprop, pred), color = "cadetblue", linewidth = 1)+
  facet_wrap(~species)+
  theme_bw()+
  xlab("Proportion industry preferred")+
  ylab("Size at 50% maturity (mm)")

# different
tan_mat <- crabpack::get_male_maturity(species = species,
                                        region = "EBS",
                                        #years = years, # set in libs/params script
                                        channel = channel)

tan_L50 <- tan_mat$model_parameters %>%
  filter(DISTRICT == "ALL") %>%
  dplyr::select(YEAR, B_EST, B_SE) %>%
  mutate(species = "Tanner")

pp <- right_join(rbind(L50, tan_L50), dat)

p1 <- pp %>% filter(species == "Snow")
# p2 <- pp %>% filter(species == "Tanner East")
# p3 <- pp %>% filter(species == "Tanner West")
p4 <- pp %>% filter(species == "Tanner") %>% mutate(survey_indpref_bio = survey_indpref_bio/1000,
                                                    df_bio = df_bio/1000)
pp2 <- rbind(p1, p4)

m1 <- gamm(B_EST ~ s(df_bio), data = na.omit(p1), correlation = corAR1(form = ~YEAR))
# m2 <- gamm(B_EST ~ s(df_bio), data = na.omit(p2), correlation = corAR1(form = ~YEAR))
# m3 <- gamm(B_EST ~ s(df_bio), data = na.omit(p3), correlation = corAR1(form = ~YEAR))
m4 <- gamm(B_EST ~ s(df_bio), data = na.omit(p4), correlation = corAR1(form = ~YEAR))


pr1 <- cbind(pred = predict(m1$gam, na.omit(p1), se = TRUE), na.omit(p1)) %>% mutate(species = "Snow")
pr4 <- cbind(pred = predict(m4$gam, na.omit(p4), se = TRUE), na.omit(p4)) %>% mutate(species = "Tanner")



ggplot()+
  geom_point(pp2, mapping = aes(df_bio, B_EST))+
  geom_ribbon(pr1 , mapping = aes(x = df_bio, ymin = pred.fit - pred.se.fit*1.96, ymax = pred.fit + pred.se.fit*1.96*1.96), fill = "cadetblue", alpha = 0.25)+
  geom_ribbon(pr4 , mapping = aes(x = df_bio, ymin = pred.fit - pred.se.fit*1.96, ymax = pred.fit + pred.se.fit*1.96*1.96), fill = "cadetblue", alpha = 0.25)+
  geom_line(pr1, mapping = aes(df_bio, pred.fit), color = "cadetblue", linewidth = 1)+
  geom_line(pr4, mapping = aes(df_bio, pred.fit), color = "cadetblue", linewidth = 1)+
  facet_wrap(~species, scales = "free")+
  theme_bw()+
  xlab("Directed fishery catch (mt)")+
  ylab("Size at 50% maturity (mm)")

# different
m1 <- gamm(B_EST ~ s(prop), data = p1, correlation = corAR1(form = ~YEAR))
# m2 <- gamm(B_EST ~ s(prop), data = p2, correlation = corAR1(form = ~YEAR))
# m3 <- gamm(B_EST ~ s(df_bio), data = p3, correlation = corAR1(form = ~YEAR))
m4 <- gamm(B_EST ~ s(prop), data = p4, correlation = corAR1(form = ~YEAR))


pr1 <- cbind(pred = predict(m1, na.omit(p1)), na.omit(p1)) %>% mutate(species = "Snow")
pr2 <- cbind(pred = predict(m2, na.omit(p2)), na.omit(p2)) %>% mutate(species = "Tanner East")



ggplot()+
  geom_point(pp, mapping = aes(prop, B_EST))+
  geom_text(pp, mapping = aes(prop, B_EST, label = YEAR))+
  geom_line(pr1, mapping = aes(prop, pred), color = "cadetblue", linewidth = 1)+
  #geom_line(pr2, mapping = aes(prop, pred), color = "cadetblue", linewidth = 1)+
  facet_wrap(~species, scales = "free")+
  theme_bw()+
  xlab("Survey industry-preferred biomass/directed fishery biomass")+
  ylab("Size at 50% maturity (mm)")

pp %>% mutate(prop = format(prop, scientific = FALSE))
