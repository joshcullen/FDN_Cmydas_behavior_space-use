
#################################################
### Compare space-use estimates among methods ###
#################################################

library(tidyverse)
library(lubridate)
library(amt)
library(move)
library(sf)
library(rnaturalearth)
library(MetBrewer)
library(units)
library(patchwork)
library(wesanderson)



#### Load the model results from each method ####

load("Processed_data/MCP_fits.RData")
load("Processed_data/KDE_fits.RData")
load("Processed_data/dBBMM_fits.RData")

dat <- read.csv('Processed_data/SSM_CRW8hr_FDN Cmydas tracks.csv')


#### Wrangle model results to match up properly ####

# Transform project for MCP to match other output
mcp <- dat.id.mcp %>%
  st_transform(crs = st_crs(contours2)) %>%
  mutate(method = 'MCP') %>%
  dplyr::select(id, level, method, geometry)

# Rename other output
kde.href <- dat.id.kde.href %>%
  mutate(method = 'KDE_href') %>%
  dplyr::select(id, level, method, geometry)
kde.hpi <- dat.id.kde.hpi %>%
  mutate(method = 'KDE_hpi') %>%
  dplyr::select(id, level, method, geometry)
dbbmm <- contours2 %>%
  rename(level = UD.level) %>%
  mutate(method = 'dBBMM') %>%
  dplyr::select(id, level, method, geometry)




#### Visually compare UDs among methods ####

brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs") %>%
  dplyr::select(-level)

ud.all <- rbind(mcp, kde.href, kde.hpi, dbbmm)


# Show all IDs, methods, and levels
ind_1_5 <- unique(dat$id)[1:5]

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat |>
              mutate(id = as.character(id)) |>
              filter(id %in% ind_1_5), aes(x, y, group = id), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = ud.all |>
            filter(id %in% ind_1_5), aes(color = method), fill = "transparent",
          linewidth = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  coord_sf(xlim = c(min(dat$x) - 50, max(dat$x) + 50),
           ylim = c(min(dat$y) - 50, max(dat$y) + 50)) +
  facet_grid(id ~ level)
#large disparity among methods for migrating IDs

# Zoom in on resident IDs
residents <- c(205542, 205544, 226069, 226071, 226072)
ggplot() +
  geom_path(data = dat |>
              filter(id %in% residents), aes(x, y, group = id), linewidth = 0.5,
            alpha = 0.15) +
  geom_sf(data = ud.all |>
            filter(id %in% residents), aes(color = method), fill = "transparent",
          linewidth = 0.5) +
  scale_color_met_d('Egypt') +
  theme_bw() +
  facet_grid(id ~ level)
#all methods relatively comparable here



############################
### Fig 5 for manuscript ###
############################

# Change labels for isopleths
ud.all$level <- ifelse(ud.all$level == 0.5, "50%", "95%")

# Plot 3 migrants and 3 residents in composite fig
mig.ids <- c(205538, 226066, 226070)
mig.plot <- ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat %>%
              filter(id %in% mig.ids), aes(x, y, group = id), linewidth = 0.5, alpha = 0.5) +
  geom_sf(data = ud.all %>%
            filter(id %in% mig.ids), aes(color = method), fill = "transparent",
          linewidth = 0.5) +
  # scale_color_manual("", values = c(MetPalettes$Juarez[[1]][c(1:3)], MetPalettes$VanGogh3[[1]][5])) +
  scale_color_met_d("Hokusai3", direction = -1) +
  labs(x = "Longitude", y = "Latitude", title = "Migratory") +
  theme_bw() +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(linewidth = 1))) +
  coord_sf(xlim = c(min(dat$x) - 100, max(dat$x) + 100),
           ylim = c(min(dat$y) - 100, max(dat$y) + 100)) +
  facet_grid(id ~ level)


res.dat <- dat %>%
  filter(id %in% c(205542, 226069, 226072)) |>
  mutate(id = as.character(id))
res.plot <- ggplot() +
  # geom_sf(data = brazil) +
  geom_path(data = res.dat, aes(x, y, group = id), linewidth = 0.5, alpha = 0.25) +
  geom_sf(data = ud.all %>%
            filter(id %in% unique(res.dat$id)), aes(color = method), fill = "transparent",
          linewidth = 0.5) +
  # scale_color_manual("", values = c(MetPalettes$Juarez[[1]][c(1:3)], MetPalettes$VanGogh3[[1]][5])) +
  scale_color_met_d("Hokusai3", direction = -1) +
  labs(x = "Longitude", y = "Latitude", title = 'Resident') +
  theme_bw() +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(linewidth = 1))) +
  coord_sf(xlim = c(min(res.dat$x) - 1, max(res.dat$x) + 1),
           ylim = c(min(res.dat$y) - 1, max(res.dat$y) + 1)) +
  facet_grid(id ~ level)



mig.plot + res.plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') +
  plot_layout(guides = 'collect') & theme(plot.tag = element_text(size = 16))

# ggsave("Figures/Fig 5.png", width = 12, height = 5, units = "in", dpi = 400)



#### Compare estimated area of space-use ####

ud.all$area <- st_area(ud.all)
ud.all$strategy <- ifelse(ud.all$id %in% residents, 'Resident', 'Migratory')


# Compare by method and UD level
iso.plot <- ggplot(ud.all, aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15,
                                                          dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3", direction = -1) +
  labs(x = "UD Isopleth") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())


# Compare by movement strategy
strategy.plot <- ggplot(ud.all %>%
         filter(level == "95%"), aes(strategy, area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15,
                                                          dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3", direction = -1) +
  labs(x = "95% UD") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

# Focus on Residents
res.area.plot <- ggplot(ud.all %>%
         filter(strategy == 'Resident'), aes(factor(level), area, color = method)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.55)) +
  geom_point(alpha = 0.7, position = position_jitterdodge(jitter.width = 0.15,
                                                          dodge.width = 0.55)) +
  scale_color_met_d("Hokusai3", direction = -1) +
  labs(x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())


############################
### Fig 6 for manuscript ###
############################

iso.plot +
  (strategy.plot + inset_element(res.area.plot, left = 0.35, right = 0.99, bottom = 0.5,
                                 top = 0.99, align_to = "panel")) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 16))

# ggsave("Figures/Fig 6.png", width = 12, height = 5, units = "in", dpi = 400)




#############################
### Fig S7 for manuscript ###
#############################

pal1 <- c(wes_palettes$Darjeeling1,
          wes_palettes$Darjeeling2,
          wes_palettes$Cavalcanti1,
          wes_palettes$Rushmore1)

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = dat %>%
              filter(!id %in% residents), aes(x, y, group = id), linewidth = 0.5,
            alpha = 0.5) +
  geom_sf(data = ud.all %>%
            filter(!id %in% residents,
                   method %in% c('KDE_hpi','dBBMM')), aes(color = id), fill = "transparent",
          linewidth = 0.5) +
  scale_color_manual(values = pal1[which(!unique(ud.all$id) %in% residents)]) +
  labs(x = "Longitude", y = "Latitude", title = "Migratory") +
  theme_bw() +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  coord_sf(xlim = c(min(dat$x) - 100, max(dat$x) + 100),
           ylim = c(min(dat$y) - 100, max(dat$y) + 100)) +
  facet_grid(method ~ level)

# ggsave("Figures/Fig S7.png", width = 10, height = 7, units = "in", dpi = 400)




#############################
### Fig S9 for manuscript ###
#############################
res.dat.all <- dat |>
  filter(id %in% residents)

ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = res.dat.all, aes(x, y, group = id), linewidth = 0.25,
            alpha = 0.5) +
  geom_sf(data = ud.all |>
            filter(id %in% residents,
                   method %in% c('KDE_hpi','dBBMM')), aes(color = id), fill = "transparent",
          linewidth = 0.5) +
  # scale_color_manual(values = pal1[which(unique(ud.all$id) %in% residents)]) +
  scale_color_viridis_d(option = "turbo") +
  labs(x = "Longitude", y = "Latitude", title = "Resident") +
  theme_bw() +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_blank()) +
  coord_sf(xlim = c(min(res.dat.all$x) - 2, max(res.dat.all$x) + 2),
           ylim = c(min(res.dat.all$y) - 2, max(res.dat.all$y) + 2)) +
  facet_grid(method ~ level)

# ggsave("Figures/Fig S9.png", width = 10, height = 7, units = "in", dpi = 400)
