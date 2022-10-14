
######################################################################################
### Fit non-parametric Bayesian movement models to track segments and observations ###
######################################################################################


library(tidyverse)
library(lubridate)
library(bayesmove)  #v0.2.1
library(sf)  #v1.0.7
library(rnaturalearth)
library(plotly)
library(furrr)
library(future)


#### Load data ####

dat_1hr <- read.csv('Processed_data/SSM_CRW1hr_FDN Cmydas tracks.csv')
dat_4hr <- read.csv('Processed_data/SSM_CRW4hr_FDN Cmydas tracks.csv')
dat_8hr <- read.csv('Processed_data/SSM_CRW8hr_FDN Cmydas tracks.csv')

glimpse(dat_1hr)

summary(dat_1hr)
summary(dat_4hr)
summary(dat_8hr)




#### Wrangle data for analysis using {bayesmove} ####

# Convert all 'date' to datetime format
dat_1hr <- dat_1hr %>%
  mutate(date = as_datetime(date))

dat_4hr <- dat_4hr %>%
  mutate(date = as_datetime(date))

dat_8hr <- dat_8hr %>%
  mutate(date = as_datetime(date))


# Calculate step lengths, turning angles, and net-squared displacement
dat_1hr_2<- prep_data(dat = dat_1hr, coord.names = c("x","y"), id = "id")
dat_4hr_2<- prep_data(dat = dat_4hr, coord.names = c("x","y"), id = "id")
dat_8hr_2<- prep_data(dat = dat_8hr, coord.names = c("x","y"), id = "id")

head(dat_8hr_2)
#since x and y are in km, steps and NSD are in km
# calculates step lengths, turning angles, net-squared displacement (NSD), and time step (dt)


# Let's double-check that all time-steps are at a single interval
table(dat_1hr_2$dt)  #yes
table(dat_4hr_2$dt)  #yes
table(dat_8hr_2$dt)  #yes


# Since we don't need to filter out obs at other time intervals, we still need to add required variables to data.frame
dat_1hr_2 <- dat_1hr_2 %>%
  group_by(id) %>%  #need to number rows separately for each ID
  mutate(time1 = 1:n(),
         obs = 1:n()) %>%
  ungroup()

dat_4hr_2 <- dat_4hr_2 %>%
  group_by(id) %>%  #need to number rows separately for each ID
  mutate(time1 = 1:n(),
         obs = 1:n()) %>%
  ungroup()

dat_8hr_2 <- dat_8hr_2 %>%
  group_by(id) %>%  #need to number rows separately for each ID
  mutate(time1 = 1:n(),
         obs = 1:n()) %>%
  ungroup()


#verify that it worked properly
dat_8hr_2 %>%
  dplyr::select(id, date, time1, obs) %>%   #select only a few cols since tibble hides time1 and obs
  split(.$id) %>%
  head()


# For direct comparison w/ HMM results, create displacement variable
dat_1hr_2$disp <- sqrt(dat_1hr_2$NSD)
dat_4hr_2$disp <- sqrt(dat_4hr_2$NSD)
dat_8hr_2$disp <- sqrt(dat_8hr_2$NSD)









#############################################
### NP Bayesian models for 1 hr time step ###
#############################################


#### Discretize data streams for models ####

# Viz density plots of each data stream
ggplot(dat_1hr_2) +
  geom_density(aes(step), fill = "cadetblue") +
  theme_bw()

ggplot(dat_1hr_2) +
  geom_density(aes(angle), fill = "firebrick") +
  theme_bw()

ggplot(dat_1hr_2) +
  geom_density(aes(disp), fill = "goldenrod") +
  theme_bw()



# Define bin limits (and number of bins)

# turning angle (naturally constrained in [0,2*pi] or [-pi,+pi])
angle.bin.lims <- c(-pi, -pi/2, -pi/4, -pi/8, -pi/16, 0, pi/16, pi/8, pi/4, pi/2, pi)  #10 bins

# step length (must be positive, but no upper bound)
step.bin.lims <- c(seq(0, 1, by = 0.2), 2, max(dat_1hr_2$step, na.rm = TRUE))  #7 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims <- c(0, 10, seq(from = 150, to = 1050, by = 150))  #8 bins


angle.bin.lims
step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc.1hr <- discrete_move_var(dat_1hr_2,
                                  lims = list(step.bin.lims, angle.bin.lims, disp.bin.lims),
                                  varIn = c("step","angle","disp"),
                                  varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc.1hr) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc.1hr) +
  geom_bar(aes(TA), fill = "firebrick") +
  theme_bw()

ggplot(dat.disc.1hr) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()






#### Fit segment-level mixed-membership model (M4) to estimate states ####

# Convert data to list by ID
dat.list.1hr <- dat.disc.1hr %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub.1hr <- purrr::map(dat.list.1hr,
                               subset,
                               select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 100000  # number of iterations for Gibbs sampler
nbins<- c(7,10,8)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg1hr<- segment_behavior(data = dat.list.sub.1hr, ngibbs = ngibbs, nbins = nbins,
                                  alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 3.5 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg1hr, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg1hr, type = "nbrks")




# Determine MAP for selecting breakpoints
MAP.est.1hr<- get_MAP(dat = dat.res.seg1hr$LML, nburn = ngibbs/2)
brkpts.1hr<- get_breakpts(dat = dat.res.seg1hr$brkpts, MAP.est = MAP.est.1hr)


# How many breakpoints estimated per ID?
apply(brkpts.1hr[,-1], 1, function(x) length(purrr::discard(x, is.na)))
brkpts.1hr




# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list.1hr, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts.1hr)

plot_breakpoints(data = dat.list.1hr, as_date = FALSE, var_names = c("SL","TA","Disp"),
                 var_labels = c("Step Length", "Turning Angle", "Displacement"),
                 brkpts = brkpts.1hr)


# Assign track segments to each ID
dat.seg.1hr<- assign_tseg(dat = dat.list.1hr, brkpts = brkpts.1hr)

head(dat.seg.1hr)








#### Cluster segments into behavioral states ####

#Select only id, tseg, and discretized data streams
dat.seg.1hr2<- dat.seg.1hr[,c("id","tseg","SL","TA","Disp")]

#Summarize observations by track segment
nbins<- c(7,10,8)
obs<- summarize_tsegs(dat = dat.seg.1hr2, nbins = nbins)
obs


set.seed(123)

# Prepare for Gibbs sampler
ngibbs<- 10000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- 9  #same as used for mixture model on observations
ndata.types<- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
dat.res.segclust.1hr<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                        ngibbs = ngibbs, nmaxclust = nmaxclust,
                                        nburn = nburn, ndata.types = ndata.types)
# takes 6 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust.1hr$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim.1hr<- extract_prop(res = dat.res.segclust.1hr, ngibbs = ngibbs, nburn = nburn,
                               nmaxclust = nmaxclust)

theta.estim_df<- theta.estim.1hr %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:all_of(nmaxclust), names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:nmaxclust

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


#Calculate mean proportions per behavior
(theta.means<- round(colMeans(theta.estim.1hr), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)  #probably 6 states



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res.seg.1hr<- get_behav_hist(dat = dat.res.segclust.1hr, nburn = nburn, ngibbs = ngibbs,
                                   nmaxclust = nmaxclust,
                                   var.names = c("Step Length","Turning Angle","Displacement"))

# Add bin lim range to each label
step.lims <- data.frame(bin.vals = cut(dat_1hr_2$step, step.bin.lims) %>%
                          levels(),
                        bin = 1:(length(step.bin.lims) - 1),
                        var = "Step Length")
angle.lims <- data.frame(bin.vals = cut(dat_1hr_2$angle, round(angle.bin.lims, 2)) %>%
                           levels(),
                         bin = 1:(length(angle.bin.lims) - 1),
                         var = "Turning Angle")
disp.lims <- data.frame(bin.vals = cut(dat_1hr_2$disp, round(disp.bin.lims, 2)) %>%
                          levels(),
                        bin = 1:(length(disp.bin.lims) - 1),
                        var = "Displacement")
lims <- rbind(step.lims, angle.lims, disp.lims)

behav.res.seg.1hr <- left_join(behav.res.seg.1hr, lims, by = c('var','bin')) %>%
  arrange(desc(var))
behav.res.seg.1hr$bin.vals <- factor(behav.res.seg.1hr$bin.vals, levels = unique(behav.res.seg.1hr$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.seg.1hr, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(7), rep("grey35", 2)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")
#actually looks like states 1-8 make sense; but states 2,3,7,8 are all foraging


# Merge states 2,3,7,8 together
tmp <- behav.res.seg.1hr %>%
  split(.$behav)
tmp[[2]]$prop <- (tmp[[2]]$prop + tmp[[3]]$prop + tmp[[7]]$prop + tmp[[8]]$prop) / 4  #calc mean of state-dependent distribs
behav.res.seg.1hr2 <- tmp[c(1:2,4:6,9)] %>%
  bind_rows()

ggplot(behav.res.seg.1hr2, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(5), rep("grey35", 2)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")


theta.estim.1hr[,2] <- rowSums(theta.estim.1hr[,c(2:3,7:8)])
theta.estim.1hr <- theta.estim.1hr[,c(1:2,4:6,9)]  #remove all other "foraging" cols


#Reformat proportion estimates for all track segments
theta.estim.long.1hr<- expand_behavior(dat = dat.seg.1hr, theta.estim = theta.estim.1hr, obs = obs, nbehav = 5,
                                       behav.names = c("Migratory","Foraging","Breeding_Encamped","Breeding_ARS","Breeding_Exploratory"),
                                       behav.order = c(3:5,2,1))

#Plot results
ggplot(theta.estim.long.1hr) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
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




#### Assign states to segments and map ####

# Convert segmented dataset into list
dat.seg.list.1hr<- df_to_list(dat = dat.seg.1hr, ind = "id")

# Merge results with original data
dat.out.1hr<- assign_behavior(dat.orig = dat.seg.1hr,
                              dat.seg.list = dat.seg.list.1hr,
                              theta.estim.long = theta.estim.long.1hr,
                              behav.names = levels(theta.estim.long.1hr$behavior))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = dat.out.1hr, aes(x=x, y=y), color="grey60", size=0.25) +
  geom_point(data = dat.out.1hr, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



# Proportion of a given state (e.g., migratory and foraging)
ggplot() +
  geom_path(data = dat.out.1hr, aes(x, y, color = Migratory, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nMigratory", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Migratory") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



ggplot() +
  geom_path(data = dat.out.1hr, aes(x, y, color = Foraging, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.1hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nForaging", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Foraging") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



dat.out.1hr %>%
  dplyr::select(-s.se) %>%  #need to remove any columns that contain only NAs
  shiny_tracks(., epsg = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")











#############################################
### NP Bayesian models for 4 hr time step ###
#############################################


#### Discretize data streams for models ####

# Viz density plots of each data stream
ggplot(dat_4hr_2) +
  geom_density(aes(step), fill = "cadetblue") +
  theme_bw()

ggplot(dat_4hr_2) +
  geom_density(aes(angle), fill = "firebrick") +
  theme_bw()

ggplot(dat_4hr_2) +
  geom_density(aes(disp), fill = "goldenrod") +
  theme_bw()



# Define bin limits (and number of bins)

# turning angle (naturally constrained in [0,2*pi] or [-pi,+pi])
angle.bin.lims <- seq(from = -pi, to = pi, by = pi/4)  #8 bins

# step length (must be positive, but no upper bound)
step.bin.lims <- c(0, 0.5, 1:5, 10, max(dat_4hr_2$step, na.rm = TRUE))  #8 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims <- seq(from = 0, to = 1050, by = 150)  #7 bins


angle.bin.lims
step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc.4hr <- discrete_move_var(dat_4hr_2,
                                  lims = list(step.bin.lims, angle.bin.lims, disp.bin.lims),
                                  varIn = c("step","angle","disp"),
                                  varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc.4hr) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc.4hr) +
  geom_bar(aes(TA), fill = "firebrick") +
  theme_bw()

ggplot(dat.disc.4hr) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()






#### Fit segment-level mixed-membership model (M4) to estimate states ####

# Convert data to list by ID
dat.list.4hr <- dat.disc.4hr %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub.4hr <- purrr::map(dat.list.4hr,
                               subset,
                               select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 100000  # number of iterations for Gibbs sampler
nbins<- c(8,8,7)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg4hr<- segment_behavior(data = dat.list.sub.4hr, ngibbs = ngibbs, nbins = nbins,
                                  alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 2 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg4hr, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg4hr, type = "nbrks")




# Determine MAP for selecting breakpoints
MAP.est.4hr<- get_MAP(dat = dat.res.seg4hr$LML, nburn = ngibbs/2)
brkpts.4hr<- get_breakpts(dat = dat.res.seg4hr$brkpts, MAP.est = MAP.est.4hr)


# How many breakpoints estimated per ID?
apply(brkpts.4hr[,-1], 1, function(x) length(purrr::discard(x, is.na)))
brkpts.4hr




# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list.4hr, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts.4hr)

plot_breakpoints(data = dat.list.4hr, as_date = FALSE, var_names = c("SL","TA","Disp"),
                 var_labels = c("Step Length", "Turning Angle", "Displacement"),
                 brkpts = brkpts.4hr)


# Assign track segments to each ID
dat.seg.4hr<- assign_tseg(dat = dat.list.4hr, brkpts = brkpts.4hr)

head(dat.seg.4hr)








#### Cluster segments into behavioral states ####

#Select only id, tseg, and discretized data streams
dat.seg.4hr2<- dat.seg.4hr[,c("id","tseg","SL","TA","Disp")]

#Summarize observations by track segment
nbins<- c(8,8,7)
obs<- summarize_tsegs(dat = dat.seg.4hr2, nbins = nbins)
obs


set.seed(123)

# Prepare for Gibbs sampler
ngibbs<- 10000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- 9  #same as used for mixture model on observations
ndata.types<- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
dat.res.segclust.4hr<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                        ngibbs = ngibbs, nmaxclust = nmaxclust,
                                        nburn = nburn, ndata.types = ndata.types)
# takes 2 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust.4hr$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim.4hr<- extract_prop(res = dat.res.segclust.4hr, ngibbs = ngibbs, nburn = nburn,
                               nmaxclust = nmaxclust)

theta.estim_df<- theta.estim.4hr %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:all_of(nmaxclust), names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:nmaxclust

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


#Calculate mean proportions per behavior
(theta.means<- round(colMeans(theta.estim.4hr), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)  #probably 5 states



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res.seg.4hr<- get_behav_hist(dat = dat.res.segclust.4hr, nburn = nburn, ngibbs = ngibbs,
                                   nmaxclust = nmaxclust,
                                   var.names = c("Step Length","Turning Angle","Displacement"))

# Add bin lim range to each label
step.lims <- data.frame(bin.vals = cut(dat_4hr_2$step, step.bin.lims) %>%
                          levels(),
                        bin = 1:(length(step.bin.lims) - 1),
                        var = "Step Length")
angle.lims <- data.frame(bin.vals = cut(dat_4hr_2$angle, round(angle.bin.lims, 2)) %>%
                           levels(),
                         bin = 1:(length(angle.bin.lims) - 1),
                         var = "Turning Angle")
disp.lims <- data.frame(bin.vals = cut(dat_4hr_2$disp, round(disp.bin.lims, 2)) %>%
                          levels(),
                        bin = 1:(length(disp.bin.lims) - 1),
                        var = "Displacement")
lims <- rbind(step.lims, angle.lims, disp.lims)

behav.res.seg.4hr <- left_join(behav.res.seg.4hr, lims, by = c('var','bin'))
behav.res.seg.4hr$bin.vals <- factor(behav.res.seg.4hr$bin.vals, levels = unique(behav.res.seg.4hr$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.seg.4hr, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(7), rep("grey35", 2)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")
#actually looks like states 1-7 make sense; but states 3-4 and 6-7 are all foraging


# Merge states 3-4,6-7 together
tmp <- behav.res.seg.4hr %>%
  split(.$behav)
tmp[[3]]$prop <- (tmp[[3]]$prop + tmp[[4]]$prop + tmp[[6]]$prop + tmp[[7]]$prop) / 4  #calc mean of state-dependent distribs
behav.res.seg.4hr2 <- tmp[c(1:3,5,8:9)] %>%
  bind_rows()

ggplot(behav.res.seg.4hr2, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 3)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")


theta.estim.4hr[,3] <- rowSums(theta.estim.4hr[,c(3:4,6:7)])
theta.estim.4hr <- theta.estim.4hr[,c(1:3,5,8:9)]  #remove all other "foraging" cols


#Reformat proportion estimates for all track segments
theta.estim.long.4hr<- expand_behavior(dat = dat.seg.4hr, theta.estim = theta.estim.4hr, obs = obs, nbehav = 4,
                                       behav.names = c("Migratory","Breeding_Encamped","Foraging","Breeding_ARS"),
                                       behav.order = c(2,4,3,1))

#Plot results
ggplot(theta.estim.long.4hr) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
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




#### Assign states to segments and map ####

# Convert segmented dataset into list
dat.seg.list.4hr<- df_to_list(dat = dat.seg.4hr, ind = "id")

# Merge results with original data
dat.out.4hr<- assign_behavior(dat.orig = dat.seg.4hr,
                              dat.seg.list = dat.seg.list.4hr,
                              theta.estim.long = theta.estim.long.4hr,
                              behav.names = levels(theta.estim.long.4hr$behavior))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = dat.out.4hr, aes(x=x, y=y), color="grey60", size=0.25) +
  geom_point(data = dat.out.4hr, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



# Proportion of a given state (e.g., migratory and foraging)
ggplot() +
  geom_path(data = dat.out.4hr, aes(x, y, color = Migratory, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nMigratory", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Migratory") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



ggplot() +
  geom_path(data = dat.out.4hr, aes(x, y, color = Foraging, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.4hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nForaging", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Foraging") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



dat.out.4hr %>%
  dplyr::select(-s.se) %>%  #need to remove any columns that contain only NAs
  shiny_tracks(., epsg = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")











#############################################
### NP Bayesian models for 8 hr time step ###
#############################################


#### Discretize data streams for models ####

# Viz density plots of each data stream
ggplot(dat_8hr_2) +
  geom_density(aes(step), fill = "cadetblue") +
  theme_bw()

ggplot(dat_8hr_2) +
  geom_density(aes(angle), fill = "firebrick") +
  theme_bw()

ggplot(dat_8hr_2) +
  geom_density(aes(disp), fill = "goldenrod") +
  theme_bw()



# Define bin limits (and number of bins)

# turning angle (naturally constrained in [0,2*pi] or [-pi,+pi])
angle.bin.lims <- seq(from = -pi, to = pi, by = pi/4)  #8 bins

# step length (must be positive, but no upper bound)
step.bin.lims <- c(seq(from = 0, to = 5, length = 6), 10, 20, max(dat_8hr_2$step, na.rm = TRUE))  #8 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims <- seq(from = 0, to = 1050, by = 150)  #7 bins


angle.bin.lims
step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc.8hr <- discrete_move_var(dat_8hr_2,
                              lims = list(step.bin.lims, angle.bin.lims, disp.bin.lims),
                              varIn = c("step","angle","disp"),
                              varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc.8hr) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc.8hr) +
  geom_bar(aes(TA), fill = "firebrick") +
  theme_bw()

ggplot(dat.disc.8hr) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()






#### Fit segment-level mixed-membership model (M4) to estimate states ####

# Convert data to list by ID
dat.list.8hr <- dat.disc.8hr %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub.8hr <- purrr::map(dat.list.8hr,
                   subset,
                   select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 100000  # number of iterations for Gibbs sampler
nbins<- c(8,8,7)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg8hr<- segment_behavior(data = dat.list.sub.8hr, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 2 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg8hr, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg8hr, type = "nbrks")




# Determine MAP for selecting breakpoints
MAP.est.8hr<- get_MAP(dat = dat.res.seg8hr$LML, nburn = ngibbs/2)
brkpts.8hr<- get_breakpts(dat = dat.res.seg8hr$brkpts, MAP.est = MAP.est.8hr)


# How many breakpoints estimated per ID?
apply(brkpts.8hr[,-1], 1, function(x) length(purrr::discard(x, is.na)))
brkpts.8hr




# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list.8hr, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts.8hr)

plot_breakpoints(data = dat.list.8hr, as_date = FALSE, var_names = c("SL","TA","Disp"),
                 var_labels = c("Step Length", "Turning Angle", "Displacement"),
                 brkpts = brkpts.8hr)


# Assign track segments to each ID
dat.seg.8hr<- assign_tseg(dat = dat.list.8hr, brkpts = brkpts.8hr)

head(dat.seg.8hr)








#### Cluster segments into behavioral states ####

#Select only id, tseg, and discretized data streams
dat.seg.8hr2<- dat.seg.8hr[,c("id","tseg","SL","TA","Disp")]

#Summarize observations by track segment
nbins<- c(8,8,7)
obs<- summarize_tsegs(dat = dat.seg.8hr2, nbins = nbins)
obs


set.seed(123)

# Prepare for Gibbs sampler
ngibbs<- 10000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- 9  #same as used for mixture model on observations
ndata.types<- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
dat.res.segclust.8hr<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                    ngibbs = ngibbs, nmaxclust = nmaxclust,
                                    nburn = nburn, ndata.types = ndata.types)
# takes 2 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust.8hr$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim.8hr<- extract_prop(res = dat.res.segclust.8hr, ngibbs = ngibbs, nburn = nburn,
                           nmaxclust = nmaxclust)

theta.estim_df<- theta.estim.8hr %>%
  as.data.frame() %>%
  pivot_longer(., cols = 1:all_of(nmaxclust), names_to = "behavior", values_to = "prop") %>%
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:nmaxclust

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(fill = "grey35", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(color = "grey35", position = position_jitter(0.1),
              alpha = 0.3) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


#Calculate mean proportions per behavior
(theta.means<- round(colMeans(theta.estim.8hr), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)  #probably 5 states



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res.seg.8hr<- get_behav_hist(dat = dat.res.segclust.8hr, nburn = nburn, ngibbs = ngibbs,
                               nmaxclust = nmaxclust,
                               var.names = c("Step Length","Turning Angle","Displacement"))

# Add bin lim range to each label
step.lims <- data.frame(bin.vals = cut(dat_8hr_2$step, step.bin.lims) %>%
                          levels(),
                        bin = 1:(length(step.bin.lims) - 1),
                        var = "Step Length")
angle.lims <- data.frame(bin.vals = cut(dat_8hr_2$angle, round(angle.bin.lims, 2)) %>%
                           levels(),
                         bin = 1:(length(angle.bin.lims) - 1),
                         var = "Turning Angle")
disp.lims <- data.frame(bin.vals = cut(dat_8hr_2$disp, round(disp.bin.lims, 2)) %>%
                          levels(),
                        bin = 1:(length(disp.bin.lims) - 1),
                        var = "Displacement")
lims <- rbind(step.lims, angle.lims, disp.lims)

behav.res.seg.8hr <- left_join(behav.res.seg.8hr, lims, by = c('var','bin'))
behav.res.seg.8hr$bin.vals <- factor(behav.res.seg.8hr$bin.vals, levels = unique(behav.res.seg.8hr$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.seg.8hr, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(7), rep("grey35", 2)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")
#actually looks like states 1-7 make sense; but states 4-7 are all foraging


# Merge states 4-7 together
tmp <- behav.res.seg.8hr %>%
  split(.$behav)
tmp[[4]]$prop <- (tmp[[4]]$prop + tmp[[5]]$prop + tmp[[6]]$prop + tmp[[7]]$prop) / 4  #calc mean of state-dependent distribs
behav.res.seg.8hr2 <- tmp[c(1:4,8:9)] %>%
  bind_rows()

ggplot(behav.res.seg.8hr2, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 3)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")


theta.estim.8hr[,4] <- rowSums(theta.estim.8hr[,4:7])


#Reformat proportion estimates for all track segments
theta.estim.long.8hr<- expand_behavior(dat = dat.seg.8hr, theta.estim = theta.estim.8hr, obs = obs, nbehav = 4,
                                   behav.names = c("Migratory", "Breeding_ARS", "Breeding_Encamped", "Foraging"),
                                   behav.order = c(3,2,4,1))

#Plot results
ggplot(theta.estim.long.8hr) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
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




#### Assign states to segments and map ####

# Convert segmented dataset into list
dat.seg.list.8hr<- df_to_list(dat = dat.seg.8hr, ind = "id")

# Merge results with original data
dat.out.8hr<- assign_behavior(dat.orig = dat.seg.8hr,
                          dat.seg.list = dat.seg.list.8hr,
                          theta.estim.long = theta.estim.long.8hr,
                          behav.names = levels(theta.estim.long.8hr$behavior))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = dat.out.8hr, aes(x=x, y=y), color="grey60", size=0.25) +
  geom_point(data = dat.out.8hr, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_fill_viridis_d("Behavior", na.value = "grey50") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



# Proportion of a given state (e.g., migratory and foraging)
ggplot() +
  geom_path(data = dat.out.8hr, aes(x, y, color = Migratory, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nMigratory", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Migratory") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



ggplot() +
  geom_path(data = dat.out.8hr, aes(x, y, color = Foraging, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out.8hr %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
  scale_color_viridis_c("Proportion\nForaging", option = "inferno", end = 0.90) +
  labs(x = "Easting", y = "Northing", title = "Foraging") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))



dat.out.8hr %>%
  dplyr::select(-s.se) %>%  #need to remove any columns that contain only NAs
  shiny_tracks(., epsg = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")





#### Export datasets for easy loading ####

save(dat.res.segclust.1hr, dat.res.segclust.4hr, dat.res.segclust.8hr,
     behav.res.seg.1hr2, behav.res.seg.4hr2, behav.res.seg.8hr2,
     theta.estim.long.1hr, theta.estim.long.4hr, theta.estim.long.8hr,
     dat.out.1hr, dat.out.4hr, dat.out.8hr,
     file = "Processed_data/bayesmove_model_fits.RData")

