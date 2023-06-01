
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


# GPS data from three individuals
dagami_43 =read.csv('data/X43.csv',sep=';', dec = ",") %>% 
  mutate( id = "43")
dagami_44 =read.csv('data/X44.csv',sep=';', dec = ",") %>% 
  mutate( id = "44")
dagami_45 =read.csv('data/X45.csv',sep=';', dec = ",") %>% 
  mutate( id = "45")

dagami = rbind(dagami_43, 
               dagami_44, 
               dagami_45) %>% 
  subset(., Date != "10/05/2023") # pre deployment test day

colony <- tibble(
  Longitude = -52.158637, 
  Latitude  = 4.644658 )%>% 
  st_as_sf(.,coords = c("Longitude","Latitude"), crs = 4326) %>% 
  st_transform(., crs = 2971)

colony_radius = 500 # units: meters

dagami  = dagami %>% 
  mutate(ts = dmy_hms(paste(dagami$Date, dagami$Time))) %>% 
  mutate(ts = ts - hms("03:00:00")) %>% 
  rename(x_ = Longitude, y_ = Latitude) %>% 
  subset(., CRC == "OK") %>% 
  dplyr::select(id, x_, y_, ts, Fix) %>% 
  st_as_sf(.,coords = c("x_","y_"), crs = 4326) %>% 
  st_transform(., crs = 2971) %>% #
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

hab = rasterize(cycle_hire_osm_projected, raster_template, field = "NIVEAU3_15",
                      fun = "max")

names(hab) = c("sol_onf")

plot(hab$sol_onf)

#   -----------------------------------------------------------------------


trk <- dagami %>% 
  make_track(x_, y_, ts,
             id = id, Fix = Fix,
             crs = 2971) %>% 
  filter(id == 45) %>% 
  tracked_from_to(from = as.POSIXct("2023-05-19 20:00"),
                  to = as.POSIXct("2024-10-01 08:00"))

trk %>% inspect()


trk = trk %>% 
  track_resample(., rate = minutes(60), tolerance = minutes(5)) %>% 
  steps_by_burst() %>% 
  random_steps(n_control = 15) %>% 
  na.omit(., cols=c("x2_","y2_")) %>% 
  st_as_sf(.,coords = c("x2_","y2_"), crs = 2971) %>% 
  mutate(dist_colony = st_distance(., colony) %>% 
           as.numeric())  %>%
  mutate(x2_ = st_coordinates(.)[, "X"],
         y2_ = st_coordinates(.)[, "Y"]) %>%
  st_drop_geometry() %>% 
  extract_covariates(hab$sol_onf) %>% 
  filter(sl_ > 0) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_)) %>% 
  filter(dist_colony > colony_radius)

trk = trk %>% 
  na.omit()




# 5. Fit a SSF an SSF, where you use elevation as the only covariate.
m1 <- trk |> fit_ssf(case_ ~ sol_onf + dist_colony  + strata(step_id_), model = TRUE)
summary(m1)


# Levels of land use
lu_levs <- levels(trk$sol_onf)

# x1
x1_lu <- data.frame(sol_onf = factor(lu_levs, levels = lu_levs),
                    dist_colony = c(mean(trk$dist_colony), mean(trk$dist_colony), 
                                    mean(trk$dist_colony), mean(trk$dist_colony)))
# x2
x2_lu <- data.frame(sol_onf = factor("Marais_interieur", levels = lu_levs),
                    dist_colony = mean(trk$dist_colony))

# log-RSS
log_rss_lu <- log_rss(m1, x1_lu, x2_lu, ci = "se")

# Plot
log_rss_lu$df |> 
  ggplot(aes(x = sol_onf_x1 , y = log_rss, ymin = lwr, ymax = upr)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) +
  xlab("Land Use Code") +
  ylab("log-RSS vs. category Marais interieur") +
  coord_cartesian(ylim = c(-13, 5)) +
  theme_bw()


