
#################################################################################
### Fit continuous-time state-space model while accounting for location error ###
#################################################################################

library(tidyverse)
library(lubridate)
library(aniMotum)
library(sf)
library(rnaturalearth)
library(terra)
library(sfarrow)
library(tictoc)
library(plotly)
library(future)
library(patchwork)
library(ggspatial)
library(wesanderson)
library(giscoR)
library(ggshadow)
library(geomtextpath)
library(ggtext)

source('R_BehavStates/helper functions.R')


#### Load data ####

dat <- read.csv('Processed_data/Cleaned_FDN Cmydas tracks.csv')

glimpse(dat); str(dat)
summary(dat)


#### Wrangle data for analysis using {aniMotum} ####

# Convert all 'Quality' values to "G" for FastGPS data and 'Date' to datetime format
dat <- dat %>%
  mutate(Quality = ifelse(Type == 'FastGPS', 'G', Quality),
         Date = as_datetime(Date))


# Rename columns for {aniMotum}
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

aniMotum::map(fit_crw_fitted,
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
  geom_path(data = dat2, aes(lon, lat, group = as.character(id)), color = 'black') +  #raw tracks
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


############################
### Fig 1 for manuscript ###
############################

# Load bathymetry
bathym <- get_elev(rast(ext(c(-43, -31, -9, -1)), crs = 'EPSG:4326',
                        res = 0.004166667),
                   maxcell = 5e8)
isobath <- terra::as.contour(bathym, levels = -200) %>%
  st_as_sf() %>%
  st_cast("LINESTRING") %>%
  slice(1) #%>%  #only keep largest isobath line
  # st_cast("POINT") %>%
  # mutate(x = st_coordinates(.)[,1],
  #        y = st_coordinates(.)[,2]) %>%
  # st_drop_geometry()

# Create color palette
pal1 <- c(wes_palettes$Darjeeling1,
          wes_palettes$Darjeeling2,
          wes_palettes$Cavalcanti1,
          wes_palettes$Rushmore1)

brazil_states<- ne_states(country = "Brazil", returnclass = 'sf')

tracks.plot <- ggplot() +
geom_sf(data = brazil_states, fill = "grey60", size = 0.3, color = "black") +
  # geom_textpath(data = isobath, aes(x = x, y = y), label = "200 m", hjust = 0.8, linewidth = 0.5) +
  geom_textsf(data = isobath, label = "200 m", hjust = 0.8, linewidth = 0.5) +
  geom_path(data = res_crw_fitted, aes(lon, lat, group = id, color = id), linewidth = 0.75,
            alpha = 0.8) +
  scale_color_manual(values = pal1, guide = "none") +
  geom_point(data = res_crw_fitted[1,], aes(lon, lat), size = 5, alpha = 0.8, shape = 21,
             fill = "gold", stroke = 1) +
  # geom_rect(aes(xmin = -32.6, xmax = -32.3, ymin = -3.95, ymax = -3.75),
  #           color = "dodgerblue3", fill = "transparent", linewidth = 1) +
  geom_text(aes(label = "b)", x = -32.1, y = -3.9), size = 4) +
  labs(x = "Longitude", y = "Latitude") +
  geom_text(aes(label = "Fernando de Noronha", x = -33, y = -3.5), size = 4, fontface = "bold") +
  geom_sf_text(data = brazil_states |>
                 filter(postal %in% c('PI','CE','RN','PB')), aes(label = name),
               fontface = "italic", size = 4, check_overlap = TRUE,
               nudge_y = c(2.5, 1, 0.25, 0)) +
  annotation_scale(location = "bl", width_hint = 0.5, style = "ticks",
                   line_col = "black", text_col = "black", line_width = 3,
                   text_cex = 1, text_face = "bold") +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  coord_sf(xlim = c(-42, -32), ylim = c(-8, -2))




# fdn.bathym <- fdn.bathym2 <- crop(bathym, ext(c(-32.6, -32.3, -3.95, -3.75)))
# fdn.bathym2[fdn.bathym2 > 0] <- NA
fdn.isobath <- st_read("Raw_data/FDN_isobaths/BATIMETRIA_SRTM_30_LINHA_SIRGAS_2000.shp") |>
  st_transform(4326) |>
  filter(PROFUNDIDA %in% c(-50,-200)) |>
  st_crop(xmin = -32.6, xmax = -32.3, ymin = -3.95, ymax = -3.75) |>
  slice(1:2)
fdn.isobath$PROFUNDIDA <- c("50 m", "200 m")
fdn.sf <- st_read_parquet('Raw_data/Brazil_land.parquet') |>
  st_transform(4326)

residents <- c(205542, 205544, 226069, 226071, 226072)

resident.plot <- ggplot() +
  # tidyterra::geom_spatraster(data = fdn.bathym2) +
  geom_sf(data = fdn.sf, fill = "grey70", size = 0.3, color = "black") +
  geom_textsf(data = fdn.isobath, aes(label = PROFUNDIDA), hjust = 0.665, linewidth = 0.5,
              text_smoothing = 0.1) +
  geom_path(data = res_crw_fitted |>
              filter(id %in% residents), aes(lon, lat, group = id, color = id), linewidth = 0.75,
            alpha = 0.8) +
  scale_color_manual(values = pal1[unique(res_crw_fitted$id) %in% residents], guide = "none") +
  geom_text(aes(label = "Fernando de Noronha", x = -32.42, y = -3.85), size = 5) +
  labs(x = "Longitude", y = "Latitude") +
  annotation_scale(location = "bl", width_hint = 0.5, style = "ticks",
                   line_col = "black", text_col = "black", line_width = 3,
                   text_cex = 1, text_face = "bold") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  coord_sf(xlim = c(-32.48, -32.4), ylim = c(-3.89, -3.84))




# background for the globe - center buffered by earth radius
ocean <- st_point(x = c(0,0)) %>%
  st_buffer(dist = 6371000) %>%
  st_sfc(crs = "+proj=ortho +lat_0=-5 +lon_0=-36 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")
ocean_df <- st_cast(ocean, "LINESTRING") %>%
  st_coordinates() %>%
  as.data.frame()

# country polygons, cut to size
world <- gisco_countries %>%
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = "+proj=ortho +lat_0=-5 +lon_0=-36 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

# get global graticule
grid <- st_graticule()
grid_crp <- st_difference(grid, st_union(gisco_countries)) %>%
  st_intersection(ocean %>% st_transform(4326)) %>% # select visible area only
  st_transform(crs = "+proj=ortho +lat_0=-5 +lon_0=-36 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")

# create bbox for tracks
bbox <- st_sfc(st_point(c(-42, -8)), st_point(c(-32, -2)), crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_transform(crs = "+proj=ortho +lat_0=-5 +lon_0=-36 +x_0=0 +y_0=0 +R=6371000 +units=m +no_defs +type=crs")


# combining several layers of shadow
# shadow.p <- ggplot() +
#   geom_glowpath(data = ocean_df,
#                 aes(X, Y, group = "L1"),
#                 shadowcolor='grey90',
#                 colour = "white",
#                 alpha = .01,
#                 shadowalpha=0.05,
#                 shadowsize = 1.8) +
#   geom_glowpath(data = ocean_df,
#                 aes(X, Y, group = "L1"),
#                 shadowcolor='grey90',
#                 colour = "white",
#                 alpha = .01,
#                 shadowalpha=0.02,
#                 shadowsize = 1) +
#   geom_glowpath(data = ocean_df,
#                 aes(X, Y, group = "L1"),
#                 shadowcolor='grey90',
#                 colour = "white",
#                 alpha = .01,
#                 shadowalpha=0.01,
#                 shadowsize = .5)

# create inset map
inset.plot <- ggplot() +
  geom_sf(data = ocean, fill = "aliceblue", color = NA) + # background first
  geom_sf(data = world, lwd = 0.1) + # now land over the oceans
  geom_sf(data = grid_crp, color = "grey80", linewidth = 0.2) +  #graticules over water
  geom_sf(data = bbox, fill = "transparent", color = "red", linewidth = 1) +
  coord_sf(expand = FALSE) +
  theme_void()



tracks.plot +
  inset_element(inset.plot, left = 0, bottom = 0.17, right = 0.45, top = 0.57,
                align_to = "plot") +
  resident.plot +
  plot_layout(nrow = 2, tag_level = "keep") +
  plot_annotation(tag_levels = list(c("a)","","b)"))) &
  theme(plot.tag = element_text(size = 22))

# ggsave("Figures/Fig 1.png", width = 9, height = 12, units = "in", dpi = 400)






############################
### Fig 2 for manuscript ###
############################

# Show time series of move persistence param (gamma) at different time steps next to maps of mapped states
# Use 2 example of migratory tracks (IDs 205540 and 41614)


res_crw_fitted <- res_crw_fitted %>%
  mutate(time.step = "Irregular")
res_crw_1hr <- res_crw_1hr %>%
  mutate(time.step = "1hr")
res_crw_4hr <- res_crw_4hr %>%
  mutate(time.step = "4hr")
res_crw_8hr <- res_crw_8hr %>%
  mutate(time.step = "8hr")

all.mods <- rbind(res_crw_fitted, res_crw_1hr, res_crw_4hr, res_crw_8hr)
all.mods$time.step <- factor(all.mods$time.step, levels = c("Irregular","1hr","4hr","8hr"))

# Create date labels for annotating both plots
date.labs <- res_crw_8hr %>%
  filter(id %in% c(205540, 41614)) |>
  mutate(date = as_date(date)) |>
  group_by(id, date) |>
  slice(1) |>
  ungroup() |>
  filter(id == 205540 & date %in% c('2021-02-20','2021-03-20','2021-04-20') |
         id == 41614 & date %in% c('2022-02-01','2022-03-01','2022-04-01')) |>
  dplyr::select(id, date, lon, lat, x, y, g) |>
  mutate(time.step = factor("Irregular", levels = levels(all.mods$time.step)))



behav.ts.plot <- ggplot() +
  geom_line(data = all.mods %>%
              filter(id == 205540 | id == 41614), aes(date, g, color = time.step),
            alpha = 0.7, linewidth = 1) +
  scale_color_manual("Time Step", values = RColorBrewer::brewer.pal(4, "Dark2")) +
  geom_rug(data = date.labs, aes(as_datetime(date)), linewidth = 1, length = unit(0.05, "npc")) +
  labs(x = 'Date', y = expression(gamma)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  facet_wrap(~id, scales = "free_x", ncol = 1)


behav.map <-
  ggplot() +
  geom_sf(data = brazil, fill = "grey60", size = 0.3, color = "black") +
  geom_path(data = all.mods %>%
              filter(id %in% c(205540, 41614)), aes(lon, lat), alpha = 0.7) +
  geom_point(data = all.mods %>%
               filter(id %in% c(205540, 41614)), aes(lon, lat, color = g)) +
  scale_color_viridis_c(expression(gamma), option = 'inferno', breaks = c(0,0.25,0.5,0.75,1),
                        limits = c(0,1)) +
  geom_textbox(data = date.labs,
                       aes(x = lon,
                           y = lat,
                           label = format(date, "%d %b %Y"),
                           hjust = case_when(id == 41614 ~ 0.5,
                                             id == 205540 & date == '2021-03-20' ~ 0.2,
                                             id == 205540 & date == '2021-02-20' ~ 0.3,
                                             TRUE ~ 0.1),
                       vjust = case_when(id == 41614 ~ 0.3,
                                         date == '2021-04-20' ~ 2.5,
                                         date == '2021-03-20' ~ 1,
                                         date == '2021-02-20' ~ 0)),
                       size = 3,
                       fill = NA,
                       box.colour = NA) +
  geom_curve(data = date.labs |>
               filter(id == 41614),
             aes(x = lon - 1.1,
                 y = lat + 0.3,
                 xend = lon,
                 yend = lat),
             curvature = -0.1,
             arrow = arrow(length = unit(0.1, "cm")),
             alpha = 0.5) +
  geom_curve(data = date.labs |>
               filter(id == 205540),
             aes(x = case_when(date == '2021-04-20' ~ lon - 0.1,
                               date == '2021-03-20' ~ lon - 0.3,
                               date == '2021-02-20' ~ lon - 0.1),
                 y = case_when(date == '2021-04-20' ~ lat - 2.1,
                               date == '2021-03-20' ~ lat - 0.3,
                               date == '2021-02-20' ~ lat + 0.3),
                 xend = lon,
                 yend = lat),
             curvature = -0.1,
             arrow = arrow(length = unit(0.1, "cm")),
             alpha = 0.5) +
    scale_x_continuous(breaks = seq(-39, -33, by = 3)) +
    scale_y_continuous(breaks = seq(-6, -4, by = 2)) +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 16)) +
  guides(color = guide_colourbar(barwidth = 15, barheight = 1)) +
  facet_grid(time.step ~ id)


behav.ts.plot + behav.map + plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 16))

# ggsave("Figures/Fig 2.png", width = 10, height = 6, units = "in", dpi = 400)





#############################
### Fig S2 for manuscript ###
#############################

# Show time series of move persistence param (gamma) at different time steps next to maps of mapped states
# Use 2 examples of resident tracks (IDs 205542 and 226072)


behav.ts.plot <- ggplot() +
  geom_line(data = all.mods %>%
              filter(id == 205542 | id == 226072), aes(date, g, color = time.step),
            alpha = 0.7, linewidth = 1) +
  scale_color_manual("Time Step", values = RColorBrewer::brewer.pal(4, "Dark2")) +
  labs(x = 'Date', y = expression(gamma)) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  facet_wrap(~id, scales = "free_x", ncol = 1)


behav.map <- ggplot() +
  # geom_sf(data = brazil, fill = "grey60", size = 0.3, color = "black") +
  geom_path(data = all.mods %>%
              filter(id %in% c(205542, 226072)), aes(lon, lat), alpha = 0.7) +
  geom_point(data = all.mods %>%
               filter(id %in% c(205542, 226072)), aes(lon, lat, color = g)) +
  scale_color_viridis_c(expression(gamma), option = 'inferno', breaks = c(0,0.25,0.5,0.75,1),
                        limits = c(0,1)) +
  # coord_sf(xlim = c(-40, -32), ylim = c(-7, -3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  scale_x_continuous(breaks = seq(-32.47, -32.41, by = 0.03), limits = c(-32.47, -32.40)) +
  scale_y_continuous(breaks = seq(-3.88, -3.86, by = 0.02), limits = c(-3.89, -3.85)) +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 16)) +
  guides(color = guide_colourbar(barwidth = 15, barheight = 1)) +
  facet_grid(time.step ~ id)


behav.ts.plot + behav.map + plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 16))

# ggsave("Figures/Fig S2.png", width = 10, height = 6, units = "in", dpi = 400)






### Export fitted tracks ###

# write.csv(res_crw_1hr, "Processed_data/SSM_CRW1hr_FDN Cmydas tracks.csv", row.names = FALSE)
# write.csv(res_crw_4hr, "Processed_data/SSM_CRW4hr_FDN Cmydas tracks.csv", row.names = FALSE)
# write.csv(res_crw_8hr, "Processed_data/SSM_CRW8hr_FDN Cmydas tracks.csv", row.names = FALSE)
# write.csv(res_crw_fitted, "Processed_data/SSM_irreg_FDN Cmydas tracks.csv", row.names = FALSE)
# save(fit_crw_fitted, fit_crw_1hr, fit_crw_4hr, fit_crw_8hr,
#      fit_crw_jmpm_fitted, fit_crw_mpm_1hr, fit_crw_mpm_4hr, fit_crw_mpm_8hr,
#      file = "Processed_data/SSM_model_fits.RData")
