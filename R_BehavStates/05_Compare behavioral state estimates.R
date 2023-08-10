
#######################################################
### Compare behavioral state estimates among models ###
#######################################################

library(tidyverse)
library(lubridate)
library(aniMotum)
library(momentuHMM)
library(bayesmove)
library(sf)
library(rnaturalearth)
library(plotly)
library(wesanderson)
library(patchwork)



#### Load the model results from each method ####

load("Processed_data/SSM_model_fits.RData")
load("Processed_data/HMM_data_and_model_fits.RData")
load("Processed_data/bayesmove_model_fits.RData")



#### Wrangle model results to compile into data.frame ####

ssm_res<- join(ssm = fit_crw_8hr,
               mpm = fit_crw_mpm_8hr,
               what.ssm = "predicted")

hmm_res <- dat_8hr_3 %>%
  mutate(state = viterbi(fit_8hr_3states))

bayes_res <- dat.out.8hr



#### Compare state-dependent distributions from HMM and M4 (bayesmove) ####

# HMM
plot(fit_8hr_3states, plotTracks = FALSE)


# M4
behav.res.seg.8hr3 <- behav.res.seg.8hr2 %>%
  mutate(behav1 = case_when(behav == 1 ~ 'Migratory',
                            behav == 2 ~ 'Breeding_ARS',
                            behav == 3 ~ 'Breeding_Encamped',
                            behav == 4 ~ 'Foraging',
                            TRUE ~ behav)) %>%
  filter(behav %in% 1:4) %>%
  mutate(across(behav1, \(x) factor(x, levels = c('Breeding_Encamped','Breeding_ARS',
                                                  'Breeding_Exploratory', 'Foraging',
                                                  'Migratory'))
  )) %>%
  mutate(across(var, \(x) factor(x, levels = unique(var))))
levels(behav.res.seg.8hr3$var) <- c("Step Length (km)", "Turning Angle (rad)",
                                    "Displacement (km)")


ggplot(behav.res.seg.8hr3, aes(x = bin.vals, y = prop, fill = behav1)) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        strip.text.y = element_text(size = 10)) +
  scale_fill_viridis_d(guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav1 ~ var, scales = "free_x")




#### Viz time series of state estimates ####

# SSM
plot(fit_crw_mpm_8hr)

# HMM
plotStates(fit_8hr_3states)

# M4
ggplot(theta.estim.long.8hr) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", linewidth = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")



#### Viz map of behavioral state estimates across methods ####

brazil <- ne_countries(scale = 10, country = "brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")



## Compare w/ focal PTT 205537

# SSM
ssm_res_205537 <- filter(ssm_res, id == 205537)

ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = ssm_res_205537, aes(x=x, y=y), color="grey60", linewidth=0.25) +
    geom_point(data = ssm_res_205537, aes(x, y, color=g), size=1.5, alpha=0.7) +
    geom_point(data = ssm_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = ssm_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_c("SSM Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(ssm_res_205537$x - 50), max(ssm_res_205537$x + 50)),
             ylim = c(min(ssm_res_205537$y - 50), max(ssm_res_205537$y + 50)))
)


# HMM
hmm_res_205537 <- hmm_res %>%
  filter(ID == 205537) %>%
  mutate(state1 = case_when(state == 1 ~ 'Breeding_Encamped',
                            state == 2 ~ 'Foraging',
                            state == 3 ~ 'Migratory'))

ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = hmm_res_205537, aes(x=x, y=y), color="grey60", linewidth=0.25) +
    geom_point(data = hmm_res_205537, aes(x, y, color=state1), size=1.5, alpha=0.7) +
    geom_point(data = hmm_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = hmm_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_d("HMM Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(hmm_res_205537$x - 50), max(hmm_res_205537$x + 50)),
             ylim = c(min(hmm_res_205537$y - 50), max(hmm_res_205537$y + 50)))
)


# M4
bayes_res_205537 <- bayes_res %>%
  filter(id == 205537)

ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = bayes_res_205537, aes(x=x, y=y), color="grey60", linewidth=0.25) +
    geom_point(data = bayes_res_205537, aes(x, y, color=behav), size=1.5, alpha=0.7) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_d("Bayesian M4 Behavior") +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(bayes_res_205537$x - 50), max(bayes_res_205537$x + 50)),
             ylim = c(min(bayes_res_205537$y - 50), max(bayes_res_205537$y + 50)))
)

ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = bayes_res_205537, aes(x=x, y=y), color="grey60", linewidth=0.25) +
    geom_point(data = bayes_res_205537, aes(x, y, color=Foraging), size=1.5, alpha=0.7) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == 1)), aes(x, y), color = "green", pch = 21,
               size = 3, stroke = 1.25) +
    geom_point(data = bayes_res_205537 %>%
                 slice(which(row_number() == n())), aes(x, y), color = "red", pch = 24,
               size = 3, stroke = 1.25) +
    scale_color_viridis_c("Prob. Foraging (M4)", option = 'inferno', end = 0.95) +
    labs(x = "Easting", y = "Northing") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    coord_sf(xlim = c(min(bayes_res_205537$x - 50), max(bayes_res_205537$x + 50)),
             ylim = c(min(bayes_res_205537$y - 50), max(bayes_res_205537$y + 50)))
)




#############################
### Fig S8 for manuscript ###
#############################

pal1 <- c(wes_palettes$Darjeeling1,
          wes_palettes$Darjeeling2,
          wes_palettes$Cavalcanti1,
          wes_palettes$Rushmore1)

m4.map.p <- ggplot() +
  geom_sf(data = brazil) +
  geom_path(data = bayes_res, aes(x, y, group = id), color = "grey60", linewidth = 0.25) +
  geom_point(data = bayes_res, aes(x, y, color = behav), size = 1.5, alpha = 0.7) +
  scale_color_viridis_d("M4 Behavioral State") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  coord_sf(xlim = c(min(bayes_res$x - 50), max(bayes_res$x + 50)),
           ylim = c(min(bayes_res$y - 50), max(bayes_res$y + 50)))


disp.ts.p <- ggplot(bayes_res) +
  geom_path(aes(date, disp, group = id, color = id)) +
  scale_color_manual("ID", values = pal1) +
  labs(x = "Date", y = "Displacement (km)") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(color = guide_legend(ncol = 2, override.aes = list(linewidth = 1)))

m4.map.p / disp.ts.p + plot_layout(heights = c(1,1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 16))

# ggsave("Figures/Fig S8.png", width = 8, height = 6, units = "in", dpi = 400)
