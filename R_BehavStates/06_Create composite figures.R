
# Create composite plots of results for manuscript

library(tidyverse)
library(sf)
library(rnaturalearth)
library(patchwork)
library(MetBrewer)
library(ggspatial)
library(patchwork)
library(geomtextpath)
library(ggtext)


#####################
### Load datasets ###
#####################

mpm.tracks <- readRDS("Processed_data/SSM_tracks.rds") |>
  mutate(model = "MPM",
         time.step = case_when(time.step == '1hr' ~ '1 hr',
                               time.step == '4hr' ~ '4 hr',
                               time.step == '8hr' ~ '8 hr',
                               TRUE ~ time.step))
hmm.tracks <- readRDS("Processed_data/HMM_tracks.rds") |>
  mutate(model = "HMM")
m4.tracks <- readRDS("Processed_data/M4_tracks.rds") |>
  mutate(model = "M4")

m4.ts <- readRDS("Processed_data/M4_time_prop.rds") |>
  mutate(model = "M4")

# Load land spatial layer
brazil <- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')



######################
### Create Figures ###
######################


### Time series plots

# Define markers for rug plot to show how state transitions compare across methods/time steps
markers <- data.frame(date = c(as_datetime("2022-02-27 08:00:00"),
                               as_datetime("2022-03-04 09:02:00"),
                               as_datetime("2022-03-14 16:00:00")),
                      marker = 1:3
                      )

# For MPM
mpm.ts.plot <-
  ggplot() +
  geom_line(data = mpm.tracks %>%
              mutate(time.step = factor(time.step,
                                        levels = c("Irregular","1 hr","4 hr","8 hr"))) %>%
              filter(id == 41614,
                     time.step %in% c("Irregular","1 hr")),
            aes(date, g, color = time.step), alpha = 0.5, linewidth = 1) +
  geom_line(data = mpm.tracks %>%
              mutate(time.step = factor(time.step,
                                        levels = c("Irregular","1 hr","4 hr","8 hr"))) %>%
              filter(model == "MPM",
                     id == 41614,
                     time.step %in% c("4 hr","8 hr")),
            aes(date, g, color = time.step), alpha = 0.8, linewidth = 1) +
  scale_color_manual("Time Step", values = RColorBrewer::brewer.pal(4, "Dark2"),
                     guide = guide_legend(override.aes = list(alpha = 1))) +
  geom_rug(data = markers, aes(date), linewidth = 1.5) +
  labs(x = 'Date', y = expression(gamma)) +
  theme_bw() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  facet_wrap(~model, scales = "free_x", ncol = 1)


# For HMM
hmm.ts.plot <-
  ggplot() +
  geom_point(data = hmm.tracks %>%
               filter(ID == 41614), aes(date, state, color = state), alpha = 0.3) +
  scale_color_manual('', values = MetPalettes$Egypt[[1]][c(1,4,3)],
                     guide = guide_legend(nrow = 2,
                                          override.aes = list(size = 3, alpha = 1))) +
  geom_rug(data = markers, aes(date), linewidth = 1.5, length = unit(0.08, "npc")) +
  labs(x = 'Date', y = 'State') +
  theme_bw() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.justification = "left",
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  facet_grid(time.step ~ model)


# For M4
m4.ts.plot <-
  ggplot(m4.ts %>%
           filter(id == 41614)) +
  geom_area(aes(x = date, y = prop, fill = behavior), color = "black", linewidth = 0.25,
            position = "fill") +
  labs(x = "Date", y = "Probability of Behavior") +
  scale_fill_manual('',
                    values = c(MetPalettes$Egypt[[1]][1:2], "black",
                               MetPalettes$Egypt[[1]][3:4]),
                    guide = guide_legend(nrow = 2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  geom_rug(data = markers, aes(date), linewidth = 1.5, length = unit(0.08, "npc")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "right",
        legend.box.margin = unit(c(0, 0, 0, 2.5), "cm")) +
  facet_grid(time.step ~ model)



### Create composite of time series plots
mpm.ts.plot / plot_spacer() / (hmm.ts.plot + m4.ts.plot) +
  plot_layout(heights = c(1, 0.1, 1.2)) +
  plot_annotation(tag_levels = "a",
                  tag_suffix = ")") &
  theme(plot.tag = element_text(size = 22))

ggsave("Figures/Fig 3.png", width = 9, height = 12, units = "in", dpi = 400)






### Maps

# Merge tracks
combined_tracks <- bind_rows(
  mpm.tracks %>%
    mutate(model = "MPM",
           time.step = case_when(time.step == '1hr' ~ '1 hr',
                                 time.step == '4hr' ~ '4 hr',
                                 time.step == '8hr' ~ '8 hr')),
  hmm.tracks %>%
    rename(id = ID) %>%
    mutate(model = "HMM"),
  m4.tracks %>%
    mutate(model = "M4")
  ) %>%
  filter(time.step != "Irregular") |>
  mutate(model = factor(model, levels = c("MPM","HMM","M4")))


# Find location closest to specified time (i.e., markers highlighted in ts plot)
markers2 <- markers |>
  expand_grid(model = c("MPM","HMM","M4"),
              time.step = c("1 hr", "4 hr", "8 hr")) |>
  rename(key_date = date) |>
  mutate(model = factor(model, levels = c("MPM","HMM","M4")))

markers.sp <- markers2 |>
  # filter(id == 41614) |>
  left_join(combined_tracks |> filter(id == 41614),
            by = join_by(model, time.step, closest(key_date >= date))) |>
  mutate(time_diff = as.duration(date - key_date)) |>
  select(key_date, date, time_diff, marker, model, time.step, lon, lat)



ggplot() +
  # 1. Base Map (Applies to all facets)
  geom_sf(data = brazil, fill = "grey60", size = 0.3, color = "black") +

  # --- 2. MPM Model (Continuous Scale: 'g') ---
  geom_point(
    data = combined_tracks %>% filter(model == "MPM",
                                      id == 41614),
    aes(lon, lat, color = g),
    alpha = 0.8, size = 0.5
  ) +
  scale_color_viridis_c(
    expression(bold("MPM"~(gamma))),
    option = 'inferno',
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    guide = guide_colorbar(title.position = "top", barwidth = 10, barheight = 1, order = 1) # Keep guide for MPM
  ) +

  # CRITICAL STEP: Reset the color scale for the next model
  ggnewscale::new_scale_color() +

  # --- 3. HMM Model (Discrete Scale: 'state') ---
  geom_point(
    data = combined_tracks %>% filter(model == "HMM",
                                      id == 41614),
    aes(lon, lat, color = state),
    alpha = 0.8, size = 0.5
  ) +
  scale_color_manual(
    'HMM',
    values = MetPalettes$Egypt[[1]][c(1, 4, 3)],
    guide = guide_legend(title.position = "top", nrow = 2, order = 2,
                         override.aes = list(size = 3)) # Define specific guide for HMM
  ) +

  # CRITICAL STEP: Reset the color scale for the next model
  ggnewscale::new_scale_color() +

  # --- 4. M4 Model (Discrete Scale: 'behav') ---
  geom_point(
    data = combined_tracks %>% filter(model == "M4",
                                      id == 41614),
    aes(lon, lat, color = behav),
    alpha = 0.8, size = 0.5
  ) +
  scale_color_manual(
    'M4',
    values = MetPalettes$Egypt[[1]], # Using different colors for differentiation
    guide = guide_legend(title.position = "top", nrow = 2, order = 3,
                         override.aes = list(size = 3)) # Define specific guide for M4
  ) +

  #--------- Add markers corresponding w/ time series plots ------------#
  geom_point(data = markers.sp, aes(lon, lat), shape = 21, size = 3, fill = "transparent") +

  # --- 5. Coordinates and Faceting ---
  scale_x_continuous(breaks = seq(-39, -33, by = 3)) +
  scale_y_continuous(breaks = seq(-6, -4, by = 2)) +
  coord_sf(xlim = c(-36, -32), ylim = c(-7, -3), expand = FALSE) +

  # Use facet_grid to create a grid of [time.step] rows by [model] columns
  facet_grid(time.step ~ model) +

  # --- 6. Theming (Unified) ---
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    # Facet strip appearance
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 10), # Added Y strip text

    # Legend control (Legends will now stack due to ggnewscale)
    legend.position = "top", # Move combined legends to bottom for clarity
    legend.box = "horizontal",
    legend.title = element_text(size = 12, face = "bold"),

    # General appearance
    panel.grid = element_blank(),

    # Aggressively reduce plot margins (use 0 for left/right for tight fit)
    plot.margin = unit(c(0.1, 0, 0, 0), "cm")
  )
ggsave("Figures/Fig 2.png", width = 8, height = 6, units = "in", dpi = 400)







# Create plot zooming in on migratory phase for MPM to better show "slow-down" phases
mpm.tracks |>
  filter(time.step == '1hr') |>
  select(-c(x, y)) |>
  rename(x = lon,
         y = lat) |>
  bayesmove::shiny_tracks(4326)



mpm.mig.rest <- mpm.tracks |>
  filter(time.step == 'Irregular',
         id %in% c(226067, 226070))

slow.labs <- data.frame(x = c(-33.2, -37),
                        y = c(-4.8, -3.8),
                        lab = c("Pelagic slow-down instances mid-migration",
                                "Slow-down during orientation toward mainland"))


ggplot() +
  geom_sf(data = brazil, fill = "grey60", size = 0.3, color = "black") +
  geom_point(
    data = mpm.mig.rest,
    aes(lon, lat, color = g),
    alpha = 0.8, size = 0.75
  ) +
  geom_path(
    data = mpm.mig.rest,
    aes(lon, lat, group = id),
    linewidth = 0
  ) +
  # Add labels per ID along tracks
  geom_textpath(data = mpm.mig.rest,
                aes(x = lon, y = lat,
                    label = paste("ID", id),
                    hjust = case_when(id == 226067 ~ 0.4,
                                      id == 226070 ~ 0.35),
                    vjust = case_when(id == 226067 ~ 1.5,
                                      id == 226070 ~ -0.25),
                    group = id),
                    text_only = TRUE) +
  # Add labels for migration slow-downs
  geom_textbox(data = slow.labs,
               aes(x = x,
                   y = y,
                   label = lab),
               size = 3,
               fill = NA,
               box.colour = NA) +
  # Add arrows for western label
  geom_curve(aes(x = -37.2,
                 y = -3.95,
                 xend = -37.457,
                 yend = -4.373 + 0.08),
             curvature = 0.1,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             alpha = 0.8) +
  geom_curve(aes(x = -37,
                 y = -3.95,
                 xend = -35.99,
                 yend = -4.715 + 0.03),
             curvature = -0.2,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             alpha = 0.8) +
  # Add arrows for eastern label
  geom_curve(aes(x = -33.5,
                 y = -4.65,
                 xend = -33.72,
                 yend = -4.337 - 0.06),
             curvature = -0.1,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             alpha = 0.8) +
  geom_curve(aes(x = -33.3,
                 y = -4.65,
                 xend = -33.609 + 0.1,
                 yend = -4 - 0.03),
             curvature = 0.2,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             alpha = 0.8) +
  scale_color_viridis_c(
    expression("Move persistence"~(gamma)~""),
    option = 'inferno',
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    guide = guide_colorbar(barwidth = 10, barheight = 1, order = 1) # Keep guide for MPM
  ) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-39, -32), ylim = c(-6, -3), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5, style = "ticks",
                   line_col = "black", text_col = "black", line_width = 3,
                   text_cex = 1, text_face = "bold") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.5), "cm"))

ggsave("Figures/Fig 4.png", width = 8, height = 4.5, units = "in", dpi = 400)
