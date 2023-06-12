

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

# GPS data 
dagami_42 = read.csv('data/X42.csv',sep=';', dec = ",") %>% 
  mutate( id = "42")
dagami_43 =read.csv('data/X43.csv',sep=';', dec = ",") %>% 
  mutate( id = "43")
dagami_44 =read.csv('data/X44.csv',sep=';', dec = ",") %>% 
  mutate( id = "44")
dagami_45 =read.csv('data/X45.csv',sep=';', dec = ",") %>% 
  mutate( id = "45")
dagami_46 =read.csv('data/X46.csv',sep=';', dec = ",") %>% 
  mutate( id = "46")
dagami_47 =read.csv('data/X47.csv',sep=';', dec = ",") %>% 
  mutate( id = "47")
dagami_48 =read.csv('data/X48.csv',sep=';', dec = ",") %>% 
  mutate( id = "48")
dagami_50 =read.csv('data/X50.csv',sep=';', dec = ",") %>% 
  mutate( id = "50")
dagami_52 =read.csv('data/X52.csv',sep=';', dec = ",") %>% 
  mutate( id = "52")
dagami_53 =read.csv('data/X53.csv',sep=';', dec = ",") %>% 
  mutate( id = "53")

dagami = rbind(dagami_42, 
               dagami_43, 
               dagami_44, 
               dagami_45, 
               dagami_46, 
               dagami_47, 
               dagami_48, 
               dagami_50, 
               dagami_52,
               dagami_53) %>% 
  subset(., Date != "10/05/2023") # journee des tests de predeploiement avec Argos


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
  subset(., NIVEAU3_15 == "Marais_interieur" | 
           NIVEAU3_15 == "Forets_inondables" )

cycle_hire_osm_projected = st_transform(onf, "EPSG:2971")
raster_template = rast(ext(cycle_hire_osm_projected), resolution = 100,
                       crs = st_crs(cycle_hire_osm_projected)$wkt)

hab = rasterize(cycle_hire_osm_projected, raster_template, field = "NIVEAU3_15",
                fun = "max")

names(hab) = c("sol_onf")


# iSSF --------------------------------------------------------------------


trk <- dagami %>% 
  make_track(x_, y_, ts,
             id = id, Fix = Fix,
             crs = 2971) %>% 
  tracked_from_to(from = as.POSIXct("2023-05-19 20:00"),
                  to = as.POSIXct("2024-10-01 08:00"))


dat1 <- trk %>% nest(data = -id)
dat2 <- dat1 %>% 
  mutate(dat.resample = map(data, ~ track_resample(., rate = minutes(60), 
                                                   tolerance = minutes(10))))%>% 
  mutate(ssf = map(dat.resample, ~ .x %>% steps_by_burst() %>% 
                     random_steps(n_control = 500) %>% 
                     extract_covariates(hab$sol_onf))) %>% 
  select(ssf, id) %>% unnest(cols = ssf) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_),
         step_id1_ = paste(id, step_id_, sep = "."))%>% 
  filter(sl_ > 0) %>% 
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_)) %>% 
  na.omit(., cols=c("x2_","y2_")) %>% 
  st_as_sf(.,coords = c("x2_","y2_"), crs = 2971) %>% 
  mutate(dist_colony = st_distance(., colony) %>% 
           as.numeric())%>% 
  mutate(x1_ = st_coordinates(.)[, "X"],
         y1_ = st_coordinates(.)[, "Y"]) %>% 
  st_drop_geometry() %>% 
  subset(dist_colony > colony_radius) %>%
  subset(dist_colony < Inf) %>%
  mutate(log_sl_ = log(sl_), sin_ta_ = sin(ta_),
         step_id1_ = paste(id, step_id_, sep = ".")) 

trk = dat2 %>% 
  na.omit()

trk$day_night = trk$t1_ %>% 
  hour() %>% 
  cut(.,
      breaks=c(-0.1, 6, 19, 24),
      labels=c('Night', 'Day','Night'))

trk_day = trk %>% 
  subset(., day_night == "Day")

trk_night = trk %>% 
  subset(., day_night == "Night")

trk = trk_night


# Levels of land use
lu_levs <- levels(trk$sol_onf)

# x1
x1_lu <- data.frame(sol_onf = factor(lu_levs, levels = lu_levs))
# x2
x2_lu <- data.frame(sol_onf = factor("Marais_interieur", levels = lu_levs))


# iSSF for each individuals -----------------------------------------------

m42 = trk %>% 
  subset(., id == 42) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r42 = m42$df %>% 
  mutate(id = as.factor(42))


m43 = trk %>% 
  subset(., id == 43) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r43 = m43$df %>% 
  mutate(id = as.factor(43))  

m44 = trk %>% 
  subset(., id == 44) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r44 = m44$df %>% 
  mutate(id = as.factor(44))  

m45 = trk %>% 
  subset(., id == 45) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r45 = m45$df %>% 
  mutate(id = as.factor(45))  

m46 = trk %>% 
  subset(., id == 46) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r46 = m46$df %>% 
  mutate(id = as.factor(46))  

m47 = trk %>% 
  subset(., id == 47) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r47 = m47$df %>% 
  mutate(id = as.factor(47))  

#m48 = trk %>% 
#  subset(., id == 48) %>% 
#  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
#  log_rss(., x1_lu, x2_lu, ci = "se") 
#r48 = m48$df %>% 
#  mutate(id = as.factor(48))  

m50 = trk %>% 
  subset(., id == 50) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r50 = m50$df %>% 
  mutate(id = as.factor(50))  

m52 = trk %>% 
  subset(., id == 52) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r52 = m52$df %>% 
  mutate(id = as.factor(52))  

m53 = trk %>% 
  subset(., id == 53) %>% 
  fit_issf(., case_ ~ sol_onf +   strata(step_id_), model = TRUE) %>% 
  log_rss(., x1_lu, x2_lu, ci = "se") 
r53 = m53$df %>% 
  mutate(id = as.factor(53))  



res = rbind(
  r42,
  r43,
  r44,
  r45,
  r46,
  r47,
  #r48,
  r50,
  r52,
  r53)


# Plot --------------------------------------------------------------------

levels(res$sol_onf_x1) = c("Forêts inondables", "Marais intérieur")

res = res %>% 
  rename(Individual = id)

res |> 
  ggplot(aes(y = sol_onf_x1 , x = log_rss,
             group = Individual, col = Individual, las=2)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lwr, xmax = upr),
                  position = position_dodge(width = 0.7), size = 0.5) +
  ylab("Habitat") +
  xlab("log-RSS versus Marais intérieur") +
  coord_cartesian(xlim = c(-1, 4.5)) +
  theme_bw()

