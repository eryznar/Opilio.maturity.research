# COMPARE GAM PREDICTIONS ----
sel <- read.csv("./Maturity research/Data/bsfrf_sel_dat.csv") %>%
  rename(SEL = selectivity, SIZE_5MM = size) %>%
  filter(year != "GAM predictions") 

sel$SEL01 <- (sel$SEL - min(sel$SEL, na.rm = TRUE)) /
  (max(sel$SEL, na.rm = TRUE) - min(sel$SEL, na.rm = TRUE))

s.gam <- gam(SEL ~ s(SIZE_5MM), data = sel, family = Gamma(link = "log"))

new.dat <- data.frame(SIZE_5MM = seq(min(sel$SIZE_5MM), max(sel$SIZE_5MM), by = 1))

pp <- cbind(new.dat, pp = predict(s.gam, new.dat, type = "response")) %>%
        mutate(type = "k=unrestricted, family=Gamma")


s.gam <- gam(SEL ~ s(SIZE_5MM, k = 5), method = "REML", data = sel, family = Gamma(link = "log"))

pp <- cbind(new.dat, predict(s.gam, new.dat, type = "response", se.fit = TRUE)) 
pp$lower <- pp$fit - pp$se.fit*1.96
pp$upper <- pp$fit + pp$se.fit*1.96



ggplot(pp, aes(SIZE_5MM, fit))+
  geom_ribbon(pp, mapping= aes(SIZE_5MM, ymin = lower, ymax = upper), fill = "cadetblue", alpha = 0.25)+
  geom_line(linewidth = 1, color = "cadetblue")+
  theme_bw()+
  ylab("Selectivity")+
  xlab("Carapace width (mm)")+
  theme(legend.position = "bottom", legend.direction = "horizontal", axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), legend.text = element_text(size=10))

ggsave("./Maturity research/Figures/gam_selectivity_preds.png", width = 7, height = 5)


# SAM and PMAT TIMESERIES ----
SAM.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   mutate(sex = "Male", type = "k=unrestricted"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   filter(YEAR >=1989) %>%
                   mutate(sex = "Female", type = "k=unrestricted"),
                 read.csv("./Maturity research/Output/SNOW_male_modeldata_k5.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   mutate(sex = "Male", type = "k=5"),
                 read.csv("./Maturity research/Output/SNOW_female_modeldata_k5.csv") %>%
                   dplyr::select(YEAR, SAM) %>%
                   mutate(sex = "Female", type = "k=5")) %>%
  filter(YEAR >=1989) %>%
  full_join(., expand.grid(YEAR = 1989:2025, sex = c("Male", "Female"), type = c("k=unrestricted", "k=5")))


ggplot(SAM.dat, aes(YEAR, SAM, color = type)) +
  geom_line() +
  geom_point() +
  ylab("Size-at-50% maturity (mm)")+
  xlab("Year")+
  #geom_smooth(method = "lm", fill = "cadetblue", color = "cadetblue", alpha = 0.2) +
  facet_wrap(~ factor(sex, levels = c("Male", "Female")),
             scales = "free_y", nrow = 1) +
  # geom_text(data = ann_df,
  #           aes(x = x, y = y, label = lab),
  #           size = 5,
  #           parse = TRUE) +
  theme_bw()+
  #ggtitle("Size-at-50% maturity")+
  theme(axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14)) -> SAM.plot

prop.dat <- rbind(read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>%
                    dplyr::select(YEAR, PROP_INDPREF) %>%
                    rename(value = PROP_INDPREF) %>%
                    mutate(cat = "Males ≥101mm", type = "k=unrestricted"),
                  read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>%
                    dplyr::select(YEAR, PMAT_5565) %>%
                    rename(value = PMAT_5565) %>%
                    mutate(cat = "Females 55-65mm", type = "k=unrestricted"),
                  read.csv("./Maturity research/Output/SNOW_female_modeldata_k5.csv") %>%
                    dplyr::select(YEAR, PMAT_5565) %>%
                    rename(value = PMAT_5565) %>% 
                    mutate(cat = "Females 55-65mm", type = "k=5"),
                  read.csv("./Maturity research/Output/SNOW_male_modeldata_k5.csv") %>%
                    dplyr::select(YEAR, PROP_INDPREF) %>%
                    rename(value = PROP_INDPREF) %>%
                    mutate(cat = "Males ≥101mm", type = "k=5")) %>%
  full_join(., expand.grid(YEAR = 1989:2025, 
                           type = unique(.$type),
                           cat= unique(.$cat)))



ggplot(prop.dat, aes(YEAR, value, color = type)) +
  geom_line() +
  geom_point() +
  ylab("Proportion mature")+
  xlab("Year")+
  #geom_smooth(method = "lm", fill = "cadetblue", color = "cadetblue", alpha = 0.2) +
  facet_wrap(~ factor(cat, levels = c("Males ≥101mm", "Females 55-65mm")),
             scales = "free_y", nrow = 1) +
  # geom_text(data = ann_df,
  #           aes(x = x, y = y, label = lab),
  #           size = 5,
  #           parse = TRUE) +
  theme_bw()+
  #ggtitle("Proportion mature")+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14)) -> pmat.plot

# Covariate comparison between selectivity curves ----
m.dat <- read.csv("./Maturity research/Output/SNOW_male_modeldata.csv") %>% mutate(type = "k=unrestricted")
m.dat.k5 <- read.csv("./Maturity research/Output/SNOW_male_modeldata_k5.csv") %>% mutate(type = "k=5")

f.dat <- read.csv("./Maturity research/Output/SNOW_female_modeldata.csv") %>% mutate(type = "k=unrestricted")
f.dat.k5 <- read.csv("./Maturity research/Output/SNOW_female_modeldata_k5.csv") %>% mutate(type = "k=5")


male.dat <- rbind(m.dat %>% dplyr::select(YEAR, LG_ABUND, COHORT_ABUND, ICE, TOCC, type), 
                  m.dat.k5 %>% dplyr::select(YEAR, LG_ABUND, COHORT_ABUND, ICE, TOCC, type)) %>%
  pivot_longer(., cols = c(2:5), names_to = "cov", values_to = "value")

ggplot(male.dat, aes(YEAR, value, color = type))+
  geom_line()+
  geom_point()+
  facet_wrap(~cov, scales = "free")+
  theme_bw()

female.dat <- rbind(f.dat %>% dplyr::select(YEAR, MALE_LG_ABUND, FEM_MAT_ABUND, FEM_COHORT_ABUND, ICE, FEM_TOCC, type), 
                    f.dat.k5 %>% dplyr::select(YEAR, MALE_LG_ABUND, FEM_MAT_ABUND, FEM_COHORT_ABUND, ICE, FEM_TOCC, type)) %>%
  pivot_longer(., cols = c(2:6), names_to = "cov", values_to = "value")

ggplot(female.dat, aes(YEAR, value, color = type))+
  geom_line()+
  geom_point()+
  facet_wrap(~cov, scales = "free")+
  theme_bw()



# Compare model covariatees and significance ----
# Male SAM
unres <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm.rda")
k5 <-  readRDS("./Maturity research/Models/SNOW_maleSAM_gamm_k5.rda")
nosel <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm_nosel.rda")
cody <- readRDS("./Maturity research/Models/SNOW_maleSAM_gamm_cody.rda")

summary(unres$gam)
summary(k5$gam)
summary(nosel$gam)
summary(cody$gam)


# Male pmat 101 exploitation
unres <- readRDS("./Maturity research/Models/SNOW_malepmat101_exploitation_gam.rda")
k5 <- readRDS("./Maturity research/Models/SNOW_malepmat101_exploitation_gam_k5.rda")
nosel <- readRDS("./Maturity research/Models/SNOW_malepmat101_exploitation_gam_nosel.rda")
cody <- readRDS("./Maturity research/Models/SNOW_malepmat101_exploitation_gam_cody.rda")

summary(unres)
summary(k5)
summary(nosel)
summary(cody)

# Female SAM
unres <- readRDS("./Maturity research/Models/SNOW_femaleSAM_gamm.rda")
k5 <- readRDS("./Maturity research/Models/SNOW_femaleSAM_gamm_k5.rda")
nosel <- readRDS("./Maturity research/Models/SNOW_femaleSAM_gamm_nosel.rda")
cody <- readRDS("./Maturity research/Models/SNOW_femaleSAM_gamm_cody.rda")

summary(unres$gam)
summary(k5$gam)
summary(nosel$gam)
summary(cody$gam)

# Female pmat5565
unres <- readRDS("./Maturity research/Models/SNOW_femalepmat5565_gam.rda")
k5 <- readRDS("./Maturity research/Models/SNOW_femalepmat5565_gam_k5.rda")
nosel <- readRDS("./Maturity research/Models/SNOW_femalepmat5565_gam_nosel.rda")
cody <- readRDS("./Maturity research/Models/SNOW_femalepmat5565_gam_cody.rda")

summary(unres)
summary(k5)
summary(nosel)
summary(cody)

