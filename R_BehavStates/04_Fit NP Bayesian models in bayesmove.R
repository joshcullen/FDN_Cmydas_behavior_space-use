
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

glimpse(dat_8hr)
summary(dat_8hr)




#### Wrangle data for analysis using {bayesmove} ####

# Convert all 'date' to datetime format
dat_8hr <- dat_8hr %>%
  mutate(date = as_datetime(date))

dat_8hr_2<- prep_data(dat = dat_8hr, coord.names = c("x","y"), id = "id")
head(dat_8hr_2)
#since x and y are in km, steps and NSD are in km
# calculates step lengths, turning angles, net-squared displacement (NSD), and time step (dt)


# Let's double-check that all time-steps are at 8 hrs (28800 s)
table(dat_8hr_2$dt)  #yes


# Since we don't need to filter out obs at other time intervals, we still need to add required variables to data.frame
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
dat_8hr_2$disp <- sqrt(dat_8hr_2$NSD)






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
step.bin.lims <- c(seq(from = 0, to = 5, length = 6), max(dat_8hr_2$step, na.rm = TRUE))  #6 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims <- seq(from = 0, to = 800, by = 200)  #4 bins


angle.bin.lims
step.bin.lims
disp.bin.lims


# Discretize data streams
dat.disc <- discrete_move_var(dat_8hr_2,
                              lims = list(step.bin.lims, angle.bin.lims, disp.bin.lims),
                              varIn = c("step","angle","disp"),
                              varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc) +
  geom_bar(aes(TA), fill = "firebrick") +
  theme_bw()

ggplot(dat.disc) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()








#### Fit observation-level mixture model to estimate states ####

# Only retain columns of discretized data streams
dat.disc.sub<- dat.disc %>%
  dplyr::select(SL, TA, Disp) %>%
  data.frame()   #cluster_obs() function crashes if trying to use 'tibble'


set.seed(2022)

# Define model params
alpha = 0.1  #prior on Dirichlet distribution
ngibbs = 10000  #number of Gibbs sampler iterations
nburn = ngibbs/2  #number of burn-in iterations
nmaxclust = 9  #number of maximum possible states (clusters) present

# Run model
dat.res.obs<- cluster_obs(dat = dat.disc.sub, alpha = alpha, ngibbs = ngibbs, nmaxclust = nmaxclust,
                      nburn = nburn)
# took 7 min to run


# Inspect traceplot of log-likelihood
plot(dat.res.obs$loglikel, type = "l")
abline(v = nburn, col = "red", lwd = 2)


## Inspect and plot results
post.seq<- (nburn + 1):ngibbs  #posterior samples

theta<- dat.res.obs$theta[post.seq,]
colnames(theta)<- 1:ncol(theta)
theta1<- colMeans(theta)
theta1
# theta1<- sort(theta1, decreasing = TRUE)
cumsum(theta1)  #possibly 3 states present; represents > 90% of all obs



# Extract bin estimates for each possible state from the `phi` matrix of the model results
behav.res.obs<- get_behav_hist(dat = dat.res.obs, nburn = nburn, ngibbs = ngibbs,
                               nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle","Displacement"))
behav.res.obs$behav<- factor(behav.res.obs$behav, levels = 1:nmaxclust)


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

behav.res.obs <- left_join(behav.res.obs, lims, by = c('var','bin'))
behav.res.obs$bin.vals <- factor(behav.res.obs$bin.vals, levels = unique(behav.res.obs$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.obs, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 5)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")
##actually looks like there's 4 states



## Assign behavioral states to observations

# Using MAP estimate, threshold of 75% assignments from posterior, and most common state
z.post<- as.data.frame(dat.res.obs$z.posterior)
z.post2<- t(apply(z.post, 1, function(x) x/sum(x)))  #calculate proportions of samples from posterior distribution assigned to each state
thresh<- 0.75  #user-defined threshold percentage for classifying a state
z.post.thresh<- apply(z.post2, 1, function(x) ifelse(max(x) > thresh, which(x > thresh), NA))
z.post.max<- apply(z.post2, 1, function(x) which.max(x))
z.map <- dat.res.obs$z.MAP

## Add states to data frame
dat.states<- dat.disc %>%
  mutate(z.map = z.map,
         z.post.thresh = z.post.thresh,
         z.post.max = z.post.max)

n.states<- 4
dat.states$z.map<- ifelse(dat.states$z.map > n.states, NA, dat.states$z.map)
dat.states$z.post.thresh<- ifelse(dat.states$z.post.thresh > n.states, NA, dat.states$z.post.thresh)
dat.states$z.post.max<- ifelse(dat.states$z.post.max > n.states, NA, dat.states$z.post.max)





# Assign names to states
dat.states2<- dat.states %>%
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                ~case_when(. == 1 ~ "Breeding_Encamped",
                           . == 2 ~ "Breeding_ARS",
                           . == 3 ~ "Foraging",
                           . == 4 ~ "Migratory",
                           is.na(.) ~ "Unclassified")
  )) %>%
  mutate(across(c('z.map','z.post.thresh','z.post.max'),
                factor, levels = c('Breeding_Encamped','Breeding_ARS','Foraging',
                                   'Migratory','Unclassified')
  ))



# Inspect number of obs assigned to each state

dat.states2 %>%   # for estimates based on MAP estimate
  group_by(z.map) %>%
  tally() %>%
  mutate(prop = n/sum(n))

dat.states2 %>%   # for estimates based on threshold on posterior
  group_by(z.post.thresh) %>%
  tally() %>%
  mutate(prop = n/sum(n))

dat.states2 %>%   # for estimates based on mode of posterior
  group_by(z.post.max) %>%
  tally() %>%
  mutate(prop = n/sum(n))



# Map results
brazil <- ne_countries(scale = 10, country = "brazil", returnclass = 'sf') %>%
  st_transform(crs = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")

# Using MAP estimates
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.map), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "MAP estimate") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)

# Using estimates w/ threshold on posterior
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.post.thresh), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "Threshold on posterior") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)


# Using estimates w/ most common state from posterior
plotly::ggplotly(
  ggplot() +
    geom_sf(data = brazil) +
    geom_path(data = dat.states2, aes(x, y, group = id), color="grey60", size=0.25) +
    geom_point(data = dat.states2, aes(x, y, fill=z.post.max), size=1.5, pch=21, alpha=0.7) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == 1)) %>%
                 ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
    geom_point(data = dat.states2 %>%
                 group_by(id) %>%
                 slice(which(row_number() == n())) %>%
                 ungroup(), aes(x, y), color = "red", pch = 24, size = 3, stroke = 1.25) +
    scale_fill_manual("Behavior",
                      values = c(viridis::viridis(4), "grey50")) +
    labs(x = "Easting", y = "Northing", title = "Most common state from posterior") +
    theme_bw() +
    theme(axis.title = element_text(size = 16),
          strip.text = element_text(size = 14, face = "bold"),
          panel.grid = element_blank()) +
    guides(fill = guide_legend(label.theme = element_text(size = 12),
                               title.theme = element_text(size = 14))) +
    coord_sf(xlim = c(min(dat.states2$x - 50), max(dat.states2$x + 50)),
             ylim = c(min(dat.states2$y - 50), max(dat.states2$y + 50)))
)


dat.states2 %>%
  mutate(across(z.map:z.post.max, as.numeric)) %>%  #can only viz time series of numeric vars
  shiny_tracks(., epsg = "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs")








#### Fit segment-level mixed-membership model to estimate states ####

# Convert data to list by ID
dat.list <- dat.disc %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub<- map(dat.list,
                   subset,
                   select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 50000  # number of iterations for Gibbs sampler
nbins<- c(6,8,4)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg1<- segment_behavior(data = dat.list.sub, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 1.5 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg1, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg1, type = "nbrks")






#### Run the segmentation model (semi-supervised via pre-specification of breakpoints) ####

#Pre-define these migratory phases
dat.list <- dat.list %>%
  map(., ~{.x %>%
      mutate(phase = case_when(disp < 6 ~ 1,
                               step > 4 ~ 2,
                               disp > 6 & step < 4 ~ 3)
      )
  })

ggplot(bind_rows(dat.list), aes(date, disp)) +
  geom_path(aes(group = id, color = factor(phase))) +
  theme_bw() +
  facet_wrap(~id, scales = "free_x")


#Find breakpoints based on 'phase'
breaks<- map(dat.list, ~find_breaks(dat = ., ind = "phase"))
breaks  #since some IDs have 0 estimated breaks and model needs at least 1 for all IDs, provide 1 fake brkpt

#All IDs need at least 1 proposed breakpoint; just create dummy location
ind <- which(lengths(breaks) == 0)
breaks[ind] <- 1
breaks


set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 50000  # number of iterations for Gibbs sampler
nbins<- c(6,8,4)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg2<- segment_behavior(data = dat.list.sub, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha, breakpt = breaks)
future::plan(future::sequential)  #return to single core
# takes 1.5 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg2, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg2, type = "nbrks")




# Determine MAP for selecting breakpoints
MAP.est1<- get_MAP(dat = dat.res.seg1$LML, nburn = ngibbs/2)
brkpts1<- get_breakpts(dat = dat.res.seg1$brkpts, MAP.est = MAP.est1)

MAP.est2<- get_MAP(dat = dat.res.seg2$LML, nburn = ngibbs/2)
brkpts2<- get_breakpts(dat = dat.res.seg2$brkpts, MAP.est = MAP.est2)

# How many breakpoints estimated per ID?
apply(brkpts1[,-1], 1, function(x) length(purrr::discard(x, is.na)))
apply(brkpts2[,-1], 1, function(x) length(purrr::discard(x, is.na)))

brkpts1
brkpts2
#looks like segmentation results were entirely unchanged by pre-specified brkpts; good sign that we reached best fit model


# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts1)

plot_breakpoints(data = dat.list, as_date = FALSE, var_names = c("SL","TA","Disp"),
                 var_labels = c("Step Length", "Turning Angle", "Displacement"),
                 brkpts = brkpts1)




# Redefine bins for SL and Disp
ggplot(dat_8hr_2) +
  geom_density(aes(step), fill = "cadetblue") +
  theme_bw()

ggplot(dat_8hr_2) +
  geom_density(aes(disp), fill = "goldenrod") +
  scale_x_continuous(breaks = seq(0, 1000, by = 150)) +
  theme_bw()


# step length (must be positive, but no upper bound)
step.bin.lims2 <- c(seq(from = 0, to = 5, length = 6), 10, max(dat_8hr_2$step, na.rm = TRUE))  #7 bins

# displacement (must be positive, but no upper bound)
disp.bin.lims2 <- seq(from = 0, to = 1050, by = 150)  #7 bins

step.bin.lims2
disp.bin.lims2


# Discretize data streams
dat.disc2 <- discrete_move_var(dat_8hr_2,
                              lims = list(step.bin.lims2, angle.bin.lims, disp.bin.lims2),
                              varIn = c("step","angle","disp"),
                              varOut = c("SL","TA","Disp"))



# Viz histograms of discretized data streams
ggplot(dat.disc2) +
  geom_bar(aes(SL), fill = "cadetblue") +
  theme_bw()

ggplot(dat.disc2) +
  geom_bar(aes(Disp), fill = "goldenrod") +
  theme_bw()


# Convert data to list by ID
dat.list2 <- dat.disc2 %>%
  df_to_list(., "id")

# Only retain id and discretized data streams
dat.list.sub2<- map(dat.list2,
                   subset,
                   select = c(id, SL, TA, Disp))



# Run the segmentation model (unsupervised)
set.seed(123)

alpha<- 1  # hyperparameter for prior (Dirichlet) distribution
ngibbs<- 50000  # number of iterations for Gibbs sampler
nbins<- c(7,8,7)  # define number of bins per data stream (in order from dat.list.sub)

progressr::handlers(progressr::handler_progress(clear = FALSE))  #to initialize progress bar
future::plan(multisession, workers = availableCores() - 2)  #run MCMC chains in parallel

dat.res.seg3<- segment_behavior(data = dat.list.sub2, ngibbs = ngibbs, nbins = nbins,
                                alpha = alpha)
future::plan(future::sequential)  #return to single core
# takes 1.5 min to run



# Trace-plots for the log marginal likelihood (LML) per ID
traceplot(data = dat.res.seg3, type = "LML")  #appears to have converged for each track

# Trace-plots for the number of breakpoints per ID
traceplot(data = dat.res.seg3, type = "nbrks")


# Determine MAP for selecting breakpoints
MAP.est3<- get_MAP(dat = dat.res.seg3$LML, nburn = ngibbs/2)
brkpts3<- get_breakpts(dat = dat.res.seg3$brkpts, MAP.est = MAP.est3)

# How many breakpoints estimated per ID?
apply(brkpts1[,-1], 1, function(x) length(purrr::discard(x, is.na)))
apply(brkpts3[,-1], 1, function(x) length(purrr::discard(x, is.na)))

brkpts1
brkpts3
#looks like segmentation results have changed a little


# Viz breakpoints w/ respect to data streams
plot_breakpoints(data = dat.list2, as_date = TRUE, var_names = c("step","angle","disp"),
                 var_labels = c("Step Length (km)", "Turning Angle (rad)", "Displacement (km)"),
                 brkpts = brkpts3)
#overall seems to do better job now

plot_breakpoints(data = dat.list2, as_date = FALSE, var_names = c("SL","TA","Disp"),
                 var_labels = c("Step Length", "Turning Angle", "Displacement"),
                 brkpts = brkpts3)



# Assign track segments to each ID
dat.seg<- assign_tseg(dat = dat.list2, brkpts = brkpts3)

head(dat.seg)








#### Cluster segments into behavioral states ####

#Select only id, tseg, and discretized data streams
dat.seg2<- dat.seg[,c("id","tseg","SL","TA","Disp")]

#Summarize observations by track segment
nbins<- c(7,8,7)
obs<- summarize_tsegs(dat = dat.seg2, nbins = nbins)
obs


set.seed(2022)

# Prepare for Gibbs sampler
ngibbs<- 10000  #number of MCMC iterations for Gibbs sampler
nburn<- ngibbs/2  #number of iterations for burn-in
nmaxclust<- 7  #same as used for mixture model on observations
ndata.types<- length(nbins)  #number of data types

# Set priors for LDA clustering model
gamma1<- 0.1
alpha<- 0.1

# Run LDA model
dat.res.segclust<- cluster_segments(dat = obs, gamma1 = gamma1, alpha = alpha,
                                    ngibbs = ngibbs, nmaxclust = nmaxclust,
                                    nburn = nburn, ndata.types = ndata.types)
# takes 2 min to run


# Check traceplot of log likelihood
plot(dat.res.segclust$loglikel, type='l', xlab = "Iteration", ylab = "Log Likelihood")
abline(v = nburn, col = "red", lwd = 2)


#Determine likely number of states from proportion assigned to each segment
theta.estim<- extract_prop(res = dat.res.segclust, ngibbs = ngibbs, nburn = nburn,
                           nmaxclust = nmaxclust)

theta.estim_df<- theta.estim %>%
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
(theta.means<- round(colMeans(theta.estim), digits = 3))

#Calculate cumulative sum
cumsum(theta.means)  #probably 5 states



# Extract bin estimates from phi matrix (for behavior distribs)
behav.res.seg<- get_behav_hist(dat = dat.res.segclust, nburn = nburn, ngibbs = ngibbs,
                               nmaxclust = nmaxclust,
                               var.names = c("Step Length","Turning Angle","Displacement"))

# Add bin lim range to each label
step.lims <- data.frame(bin.vals = cut(dat_8hr_2$step, step.bin.lims2) %>%
                          levels(),
                        bin = 1:(length(step.bin.lims2) - 1),
                        var = "Step Length")
angle.lims <- data.frame(bin.vals = cut(dat_8hr_2$angle, round(angle.bin.lims, 2)) %>%
                           levels(),
                         bin = 1:(length(angle.bin.lims) - 1),
                         var = "Turning Angle")
disp.lims <- data.frame(bin.vals = cut(dat_8hr_2$disp, round(disp.bin.lims2, 2)) %>%
                          levels(),
                        bin = 1:(length(disp.bin.lims2) - 1),
                        var = "Displacement")
lims <- rbind(step.lims, angle.lims, disp.lims)

behav.res.seg <- left_join(behav.res.seg, lims, by = c('var','bin'))
behav.res.seg$bin.vals <- factor(behav.res.seg$bin.vals, levels = unique(behav.res.seg$bin.vals))

# Plot state-dependent distributions
ggplot(behav.res.seg, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(6), rep("grey35", 1)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")
#actually looks like states 1-6 make sense; but states 3, 5, and 6 are all foraging


# Merge states 3, 5, and 6 together
tmp <- behav.res.seg %>%
  split(.$behav)
tmp[[3]]$prop <- (tmp[[3]]$prop + tmp[[5]]$prop + tmp[[6]]$prop) / 3  #calc mean of state-dependent distribs
behav.res.seg2 <- tmp[c(1:4,7)] %>%
  bind_rows()

ggplot(behav.res.seg2, aes(x = bin.vals, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12, angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = c(viridis::viridis(4), rep("grey35", 1)), guide = 'none') +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  facet_grid(behav ~ var, scales = "free_x")


theta.estim[,3] <- theta.estim[,3] + theta.estim[,5] + theta.estim[,6]


#Reformat proportion estimates for all track segments
theta.estim.long<- expand_behavior(dat = dat.seg, theta.estim = theta.estim, obs = obs, nbehav = 4,
                                   behav.names = c("Migratory", "Breeding_Encamped", "Foraging",
                                                   "Breeding_ARS"),
                                   behav.order = c(2,4,3,1))

#Plot results
ggplot(theta.estim.long) +
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
dat.seg.list<- df_to_list(dat = dat.seg, ind = "id")

# Merge results with original data
dat.out<- assign_behavior(dat.orig = dat.seg,
                          dat.seg.list = dat.seg.list,
                          theta.estim.long = theta.estim.long,
                          behav.names = levels(theta.estim.long$behavior))

# Map dominant behavior for all IDs
ggplot() +
  geom_path(data = dat.out, aes(x=x, y=y), color="grey60", size=0.25) +
  geom_point(data = dat.out, aes(x, y, fill=behav), size=1.5, pch=21, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out %>%
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
  geom_path(data = dat.out, aes(x, y, color = Migratory, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out %>%
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
  geom_path(data = dat.out, aes(x, y, color = Foraging, group = id), size=0.5, alpha=0.7) +
  geom_point(data = dat.out %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3, stroke = 1.25) +
  geom_point(data = dat.out %>%
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





#### Export datasets for easy loading ####

save(behav.res.seg2, theta.estim.long, dat.out, dat.res.seg3, dat.res.segclust,
     file = "Processed_data/bayesmove_model_fits.RData")
