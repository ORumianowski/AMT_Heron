
#----------- Module 09 -- Multiple animals ------------X

# Load packages ----
library(tidyverse)
library(mvtnorm)
library(amt)
library(purrr)
library(lubridate)
library(sf)
library(terra)
library(glmmTMB)
library(broom)


# 1. Load data ----

dagami_43 =read.csv('data/X43.csv',sep=';', dec = ",") %>% 
  mutate( id = "43")
dagami_44 =read.csv('data/X44.csv',sep=';', dec = ",") %>% 
  mutate( id = "44")
dagami_45 =read.csv('data/X45.csv',sep=';', dec = ",") %>% 
  mutate( id = "45")

dagami = rbind(dagami_43, 
               dagami_44, 
               dagami_45) %>% 
  subset(., Date != "10/05/2023") # pre" deployment test day

colony <- tibble(
  Longitude = -52.158637, 
  Latitude  = 4.644658 )%>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 4326) %>% 
  st_transform(., crs = 2971)

rayon_colony = units::set_units(500, m)

dagami  = dagami %>% 
  mutate(ts = dmy_hms(paste(dagami$Date, dagami$Time))) %>% 
  mutate(ts = ts - hms("03:00:00")) %>%  # mise a l heure de Cayenne
  rename(x_ = Longitude, y_ = Latitude) %>% 
  subset(., CRC == "OK") %>% 
  dplyr::select(id, x_, y_, ts, Fix) %>% 
  st_as_sf(.,coords = c("x_","y_"), crs = 4326) %>% 
  st_transform(., crs = 2971) %>% 
  mutate(x_ = st_coordinates(.)[, "X"],
         y_ = st_coordinates(.)[, "Y"]) %>%
  st_drop_geometry()



# Habitat data 

onf = st_read("ocs_onf_2015/oc_sol_2015.shp") %>% 
  dplyr::select(NIVEAU3_15) %>% 
  mutate(NIVEAU3_15 = fct_collapse(NIVEAU3_15,
                                   Mangroves = c("318"),
                                   Marais_interieur = c("411"),
                                   Forets_inondables = c("317"), 
                                   Forets_terrestres = c("3151", "3152", "3161", "3162", "3154"))) %>% 
  subset(., NIVEAU3_15 == "Mangroves" | NIVEAU3_15 == "Marais_interieur" | 
           NIVEAU3_15 == "Forets_inondables" | NIVEAU3_15 == "Forets_terrestres")

cycle_hire_osm = onf
cycle_hire_osm_projected = st_transform(cycle_hire_osm, "EPSG:2971")
raster_template = rast(ext(cycle_hire_osm_projected), resolution = 100,
                       crs = st_crs(cycle_hire_osm_projected)$wkt)

hab_terra = rasterize(cycle_hire_osm_projected, raster_template, field = "NIVEAU3_15",
                      fun = "max")

names(hab_terra) = c("sol_onf")



# iSSF --------------------------------------------------------------------

trk <- dagami %>% 
  make_track(x_, y_, ts,
             id = id, Fix = Fix,
             crs = 2971) %>% 
  tracked_from_to(from = as.POSIXct("2023-05-19 20:00"),
                  to = as.POSIXct("2024-10-01 08:00"))


dat1 <- trk %>% nest(data = -id)

dat2 <- dat1 %>% 
  mutate(dat.resample = map(data, ~ track_resample(., rate = minutes(60), tolerance = minutes(10))))

dat.ssf <- dat2 %>% 
  mutate(ssf = map(dat.resample, ~ .x %>% steps_by_burst() %>% 
                     random_steps(n_control = 10) %>% extract_covariates(hab_terra$sol_onf))) 

dat.ssf2 <- dat.ssf %>% select(ssf, id) %>% unnest(cols = ssf) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_),
         step_id1_ = paste(id, step_id_, sep = "."))


dat.ssf3 = dat.ssf2 %>% 
  filter(sl_ > 0) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_)) %>% 
  na.omit(., cols=c("x2_","y2_")) %>% 
  st_as_sf(.,coords = c("x2_","y2_"), crs = 2971) %>% 
  mutate(dist_colony = st_distance(., colony))%>% 
  mutate(x1_ = st_coordinates(.)[, "X"],
         y1_ = st_coordinates(.)[, "Y"]) %>% 
  st_drop_geometry() %>% 
  subset(dist_colony > rayon_colony) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_),
         step_id1_ = paste(id, step_id_, sep = ".")) 

dat.ssf3 = dat.ssf3 %>% 
  na.omit()


# First simple model --------------------------------------------------------

m.ssf0 <- fit_ssf(case_ ~ 
                    dist_colony + sol_onf +
                    + sl_ + log_sl_ +
                    strata(step_id1_), 
                  
                  data = dat.ssf3)

m.ssf1 <- glmmTMB(case_ ~ -1 
                  + dist_colony + sol_onf  
                    + sl_ + log_sl_ + 
                    strata(step_id1_),
                  
                  family = poisson(), 
                  data = dat.ssf3, doFit = FALSE)

m.ssf1$parameters$theta[1] <- log(1e3)
m.ssf1$mapArg <- list(theta = factor(c(NA)))
m.ssf1 <- glmmTMB:::fitTMB(m.ssf1)

coef(m.ssf0)
fixef(m.ssf1) #not the same value...

confint(m.ssf0$model)
confint(m.ssf1)

vcov(m.ssf0$model)
vcov(m.ssf1)

summary(m.ssf0) 
summary(m.ssf1)


# Second simple model --------------------------------------------------------

m.ssf0 <- fit_ssf(case_ ~ dist_colony + sol_onf +
                    + sl_ + log_sl_ +
                    strata(step_id1_), data = dat.ssf3)

m.ssf1 <- glmmTMB(case_ ~ -1 + dist_colony + sol_onf  
                  + sl_ + log_sl_ + 
                    strata(step_id1_) +
                    (1 | step_id1_) + (0 + dist_colony|id),
                  family = poisson(), 
                  data = dat.ssf3, doFit = FALSE)

m.ssf1$parameters$theta[1] <- log(1e3)
m.ssf1$mapArg <- list(theta = factor(c(NA, 1)))
m.ssf1 <- glmmTMB:::fitTMB(m.ssf1)

coef(m.ssf0)
fixef(m.ssf1)

confint(m.ssf0$model)
confint(m.ssf1)

vcov(m.ssf0$model)
vcov(m.ssf1)

summary(m.ssf0) 
summary(m.ssf1)

# Complete model ----------------------------------------------------------

m.ssf0 <- fit_ssf(case_ ~ dist_colony + sol_onf +
                    + sl_ + log_sl_ +
                    strata(step_id1_), data = dat.ssf3)

m.ssf1 <- glmmTMB(case_ ~ -1 + dist_colony + sol_onf  
                  + sl_ + log_sl_ + 
                    strata(step_id1_) +
                    (1 | step_id1_) + 
                    (0 + dist_colony|id) + (0 + sol_onf | id) ,
                  family = poisson(), 
                  data = dat.ssf3, doFit = FALSE)

m.ssf1$parameters$theta[1] <- log(1e3)
m.ssf1$mapArg <- list(theta = factor(c(NA, 1:2)))
m.ssf1 <- glmmTMB:::fitTMB(m.ssf1) # Error : from the structure of mapArg ?

coef(m.ssf0)
fixef(m.ssf1)

confint(m.ssf0$model)
confint(m.ssf1)

vcov(m.ssf0$model)
vcov(m.ssf1)

summary(m.ssf0) 
summary(m.ssf1)




