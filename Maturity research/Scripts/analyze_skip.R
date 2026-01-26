# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Process specimen
sh3.chela<-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION ==3, SEX == 1, is.na(CHELA_HEIGHT) == FALSE) %>%
  filter(SPECIES == "SNOW", SIZE >= 35 & SIZE <= 135) %>%
  mutate(LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT),
         CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  #filter(N_chela > 40) %>% #necessary for consecutive size correlations to run
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, MATURE) %>%
  filter(YEAR != 2012)


# plot
ggplot()+
  geom_bar(sh3.chela, mapping = aes(YEAR), stat = "count") +
  theme_bw()

ggplot()+
  geom_bar(sh3.chela, mapping = aes(SIZE_5MM), stat = "count") +
  theme_bw()

ggplot()+
  geom_bar(sh3.chela, mapping = aes(SIZE_5MM), stat = "count") +
  theme_bw()+
  facet_wrap(~YEAR, scales = "free_y")+
  ggtitle("SH3 males")


# calculate proportion skip molt
skip <- sh3.chela %>%
  group_by(YEAR) %>%
  reframe(
    n_total = n(),
    n_imm   = sum(MATURE == 0),
    p_imm   = n_imm / n_total,
    se_imm  = sqrt(p_imm * (1 - p_imm) / n_total),
    lcl_95  = p_imm - 1.96 * se_imm,
    ucl_95  = p_imm + 1.96 * se_imm
  ) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR), by = 1)))


ggplot(skip, aes(YEAR, p_imm)) +
  geom_point() +
  geom_errorbar(aes(ymin = lcl_95, ymax = ucl_95), width = 0.2) +
  ylab("Proportion skip-molt") +
  xlab("Year")+
  theme_bw()+
  geom_line()
