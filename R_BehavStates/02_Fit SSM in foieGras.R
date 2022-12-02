
#################################################################################
### Fit continuous-time state-space model while accounting for location error ###
#################################################################################

library(tidyverse)
library(lubridate)
library(foieGras)  #v1.0-7
library(sf)  #v1.0.7
library(rnaturalearth)
library(tictoc)
library(plotly)
library(future)
library(patchwork)


#### Load data ####

dat <- read.csv('Processed_data/Cleaned_FDN Cmydas tracks.csv')

glimpse(dat); str(dat)
summary(dat)


#### Wrangle data for analysis using {foieGras} ####

# Convert all 'Quality' values to "G" for FastGPS data and 'Date' to datetime format
dat <- dat %>%
  mutate(Quality = ifelse(Type == 'FastGPS', 'G', Quality),
         Date = as_datetime(Date))


# Rename columns for {foieGras}
dat2<- dat %>%
  rename(id = Ptt, date = Date, lc = Quality, lon = Longitude, lat = Latitude,
         eor = Error.Ellipse.orientation, smaj = Error.Semi.major.axis,
         smin = Error.Semi.minor.axis) %>%
  dplyr::select(id, date, lc, lon, lat, smaj, smin, eor)  #reorders and subsets the columns

glimpse(dat2)



#### Inspect time steps of transmissions for making predictions ####

tmp <- dat2 %>%
  split(.$id) %>%
  purrr::map(., ~mutate(.x,
                        dt = difftime(c(date[-1], NA),
                                      date,
                                      units = "secs") %>%
                          as.numeric())
             ) %>%
  bind_rows()


ggplot(tmp, aes(date, dt)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")
# it looks like there are some very large gaps in the tracks for PTTs 41587, 41588; let's remove these points that occur before (or after) depending on the ID


# Find rows where 1st dt of an ID is > 7 d (604800 secs)
ind <- tmp %>%
  group_by(id) %>%
  filter(dt > 3600*24*7) %>%  #what's the first observation w/ time step longer than 1 week?
  slice(1)

# Remove obs before or after a gap of 7 days
dat3 <- dat2 %>%
  filter(!(id == 41587 & date <= ind$date[1])) %>%  #for 41587, keep all points after this datetime
  filter(!(id == 41588 & date >= ind$date[2]))  #for41588, keep all points before this datetime



# Determine primary time step
ggplot(tmp) +
  geom_histogram(aes(dt), binwidth = 3600) +
  theme_bw() +
  xlim(0,3600*24)

table(tmp$dt) %>%
  sort(decreasing = TRUE)

tmp %>%
  filter(!(id == 41587 & date <= ind$date[1])) %>%
  filter(!(id == 41588 & date >= ind$date[2])) %>%
  group_by(id) %>%
  summarize(mean = mean(dt, na.rm = TRUE),
            median = median(dt, na.rm = TRUE))
# Mean/median time step is ~1 hr




#################
#### Run SSM ####
#################

# Change `id` to character to avoid problems during model runs
dat3$id <- as.character(dat3$id)



#### Account for location error at observed irregular time interval ####

# Estimate 'true' locations on irregular sampling interval (by setting `time.step = NA`)
tic()
fit_crw_fitted <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = NA,
                          control = ssm_control(verbose = 1))
toc()  #took 1 min where time.step = NA

print(fit_crw_fitted)  #all indiv. models converged


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_fitted, what = "fitted", type = 1, ask = TRUE)
plot(fit_crw_fitted, what = "fitted", type = 2, ask = TRUE)

foieGras::map(fit_crw_fitted,
    what = "fitted",
    by.id = TRUE)


# Estimate behavioral state (i.e., move persistence; gamma)
# Joint move persistence model ('jmpm') uses hierarchical approach across IDs
tic()
fit_crw_jmpm_fitted <- fit_mpm(fit_crw_fitted, what = "fitted", model = "jmpm",
                              control = mpm_control(verbose = 1))
toc()  #took 5 min to fit

print(fit_crw_jmpm_fitted)  #model converged
plot(fit_crw_jmpm_fitted)


# Grab results and plot
res_crw_fitted<- join(ssm = fit_crw_fitted,
                      mpm = fit_crw_jmpm_fitted,
                      what.ssm = "fitted")


# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat2, aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res_crw_fitted, aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -1))


# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_fitted, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_fitted, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)


#interactively explore using projected coordinates (World Mercator Projection; EPSG:3395, units = km)
res_crw_fitted %>%
  dplyr::select(id, date, x, y, s, g) %>%
  data.frame() %>%  #current version of shiny_tracks won't work w/ tibble format or when a column has all NAs
  bayesmove::shiny_tracks(., "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")


#interactively explore using lat/long (EPSG:4326)
res_crw_fitted %>%
  dplyr::select(-c(x, y, s.se)) %>%
  rename(x = lon, y = lat) %>%
  dplyr::select(id, date, x, y, s, g) %>%
  data.frame() %>%  #current version of shiny_tracks won't work w/ tibble format or when a column has all NAs
  bayesmove::shiny_tracks(., 4326)


# Check model fit w/ diagnostic "one-step-ahead residuals"; takes long time!
tic()
gof_irreg <- osar(fit_crw_fitted)
toc()  # takes ~1.7 hrs to run


par(ask = TRUE)
for (k in unique(gof_irreg$id)) {
  tmp <- gof_irreg %>%
    filter(id == k)
  print((plot(tmp, type = "ts") | plot(tmp, type = "qq")) /
          (plot(tmp, type = "acf") | plot_spacer()))
}
par(ask = FALSE)




#### Account for location error at regularized time interval with correlated random walk ####

### Estimate 'true' locations on regular sampling interval of 1 hr ###
tic()
fit_crw_1hr <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = 1,
                       control = ssm_control(verbose = 1))
toc()  #took 2 min to fit model

print(fit_crw_1hr)  #all models converged


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_1hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw_1hr, what = "predicted", type = 2, ask = TRUE)



# Estimate behavioral state (i.e., move persistence; gamma)
# Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
fit_crw_mpm_1hr <- fit_mpm(fit_crw_1hr, what = "predicted", model = "mpm",
                           control = mpm_control(verbose = 1))
toc()  #took 3 min to fit
#if issues w/ model not converging, try changing time step of SSM and re-running

print(fit_crw_mpm_1hr)  #all models converged
plot(fit_crw_mpm_1hr)


# Check GOF
tic()
gof_1hr <- osar(fit_crw_1hr)
toc()  # takes 2.5 hrs to run


par(ask = TRUE)
for (k in unique(gof_1hr$id)) {
  tmp <- gof_1hr %>%
    filter(id == k)
  print((plot(tmp, type = "ts") | plot(tmp, type = "qq")) /
          (plot(tmp, type = "acf") | plot_spacer()))
}
par(ask = FALSE)





### Estimate 'true' locations on regular sampling interval of 4 hrs ###
tic()
fit_crw_4hr <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = 4,
                       control = ssm_control(verbose = 1))
toc()  #took 73 sec to fit model

print(fit_crw_4hr)  #all models converged


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_4hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw_4hr, what = "predicted", type = 2, ask = TRUE)



# Estimate behavioral state (i.e., move persistence; gamma)
# Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
fit_crw_mpm_4hr <- fit_mpm(fit_crw_4hr, what = "predicted", model = "mpm",
                           control = mpm_control(verbose = 1))
toc()  #took 42 sec to fit
#if issues w/ model not converging, try changing time step of SSM and re-running

print(fit_crw_mpm_4hr)  #all models converged
plot(fit_crw_mpm_4hr)



# Check GOF
tic()
gof_4hr <- osar(fit_crw_4hr)
toc()  # takes 1.9 hrs to run


par(ask = TRUE)
for (k in unique(gof_4hr$id)) {
  tmp <- gof_4hr %>%
    filter(id == k)
  print((plot(tmp, type = "ts") | plot(tmp, type = "qq")) /
          (plot(tmp, type = "acf") | plot_spacer()))
}
par(ask = FALSE)





### Estimate 'true' locations on regular sampling interval of 8 hrs ###
tic()
fit_crw_8hr <- fit_ssm(dat3, vmax = 3, model = "crw", time.step = 8,
                          control = ssm_control(verbose = 1))
toc()  #took 1 min to fit model

print(fit_crw_8hr)  #all models converged


# Viz the filtered outliers (gold), raw observations (blue), and estimated locations (red), along with associated uncertainty (red shading)
plot(fit_crw_8hr, what = "predicted", type = 1, ask = TRUE)
plot(fit_crw_8hr, what = "predicted", type = 2, ask = TRUE)



# Estimate behavioral state (i.e., move persistence; gamma)
# Individual move persistence model ('mpm') estimates behavioral states separately across IDs
tic()
fit_crw_mpm_8hr <- fit_mpm(fit_crw_8hr, what = "predicted", model = "mpm",
                              control = mpm_control(verbose = 1))
toc()  #took 35 sec to fit
#if issues w/ model not converging, try changing time step of SSM and re-running

print(fit_crw_mpm_8hr)
plot(fit_crw_mpm_8hr)

# Check GOF
tic()
gof_8hr <- osar(fit_crw_8hr)
toc()  # takes ~1.8 hrs to run


par(ask = TRUE)
for (k in unique(gof_8hr$id)) {
  tmp <- gof_8hr %>%
    filter(id == k)
  print((plot(tmp, type = "ts") | plot(tmp, type = "qq")) /
    (plot(tmp, type = "acf") | plot_spacer()))
}
par(ask = FALSE)




#### Grab results and plot ####
res_crw_1hr<- join(ssm = fit_crw_1hr,
                   mpm = fit_crw_mpm_1hr,
                   what.ssm = "predicted")

res_crw_4hr<- join(ssm = fit_crw_4hr,
                   mpm = fit_crw_mpm_4hr,
                   what.ssm = "predicted")

res_crw_8hr<- join(ssm = fit_crw_8hr,
                      mpm = fit_crw_mpm_8hr,
                      what.ssm = "predicted")


# Compare raw tracks vs fitted tracks (for adults tagged at Fernando de Noronha)
ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat2, aes(lon, lat, group = id), color = 'black') +  #raw tracks
  geom_path(data = res_crw_1hr, aes(lon, lat, group = id), color = "blue") +  #modeled tracks
  geom_path(data = res_crw_4hr, aes(lon, lat, group = id), color = "red") +  #modeled tracks
  geom_path(data = res_crw_8hr, aes(lon, lat, group = id), color = "green") +  #modeled tracks
  theme_bw() +
  facet_wrap(~id) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -1))

# Viz modeled tracks together
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_8hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)



# Viz modeled tracks together w/ behavior plotted (1 hr)
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_1hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted (4 hr)
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_4hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)

# Viz modeled tracks together w/ behavior plotted (8 hr)
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_point(data = res_crw_8hr, aes(lon, lat, group = id, color = g), size = 0.75, alpha = 0.8) +
    scale_color_viridis_c(option = "inferno") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)


# Compare predicted tracks against tracks fitted at irregular time step
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil, fill = "grey60") +
    geom_path(data = res_crw_fitted, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.3) +
    geom_path(data = res_crw_8hr, aes(lon, lat, group = id, color = id), size = 0.75, alpha = 0.8) +
    scale_color_viridis_d() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))
)









### Compare estimates across the different approaches ###

# Behavioral states
ggplot() +
  geom_line(data = res_crw_fitted, aes(date, g, color = "CRW_Irregular"), alpha = 0.7) +
  geom_line(data = res_crw_1hr, aes(date, g, color = "CRW_1hr"), alpha = 0.7) +
  geom_line(data = res_crw_4hr, aes(date, g, color = "CRW_4hr"), alpha = 0.7) +
  geom_line(data = res_crw_8hr, aes(date, g, color = "CRW_8hr"), alpha = 0.7) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")

## Of these 4 different results, the CRW model w/ 8 hr prediction time step seems to give best broad scale state estimates; ecological interpretability at coarse-scale declines w/ decreasing time step


# Track relocations
ggplot() +
  geom_path(data = res_crw_fitted, aes(lon, lat, color = "CRW_Irregular"), alpha = 0.7) +
  geom_path(data = res_crw_1hr, aes(lon, lat, color = "CRW_1hr"), alpha = 0.7) +
  geom_path(data = res_crw_4hr, aes(lon, lat, color = "CRW_4hr"), alpha = 0.7) +
  geom_path(data = res_crw_8hr, aes(lon, lat, color = "CRW_8hr"), alpha = 0.7) +
  scale_color_manual(values = RColorBrewer::brewer.pal(4, "Dark2")) +
  theme_bw() +
  facet_wrap(~id, scales = "free")

## Estimated tracks are comparable across all time steps






### Export fitted tracks ###

write.csv(res_crw_1hr, "Processed_data/SSM_CRW1hr_FDN Cmydas tracks.csv", row.names = FALSE)
write.csv(res_crw_4hr, "Processed_data/SSM_CRW4hr_FDN Cmydas tracks.csv", row.names = FALSE)
write.csv(res_crw_8hr, "Processed_data/SSM_CRW8hr_FDN Cmydas tracks.csv", row.names = FALSE)
write.csv(res_crw_fitted, "Processed_data/SSM_irreg_FDN Cmydas tracks.csv", row.names = FALSE)
