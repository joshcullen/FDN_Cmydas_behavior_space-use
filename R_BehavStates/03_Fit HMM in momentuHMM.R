
#############################################
### Fit discrete-time hidden Markov model ###
#############################################

library(tidyverse)
library(lubridate)
library(momentuHMM)  #v1.5.4
library(sf)  #v1.0.7
library(tictoc)
library(plotly)
library(rnaturalearth)
library(MetBrewer)
library(patchwork)

source("R_BehavStates/helper functions.R")


#### Load data ####

dat_1hr <- read.csv('Processed_data/SSM_CRW1hr_FDN Cmydas tracks.csv')
dat_4hr <- read.csv('Processed_data/SSM_CRW4hr_FDN Cmydas tracks.csv')
dat_8hr <- read.csv('Processed_data/SSM_CRW8hr_FDN Cmydas tracks.csv')

glimpse(dat_1hr)

summary(dat_1hr)
summary(dat_4hr)
summary(dat_8hr)



#### Wrangle data for analysis using {momentuHMM} ####

# Convert all 'Date' to datetime format and change name of 'id' column
dat_1hr <- dat_1hr %>%
  mutate(date = as_datetime(date)) %>%
  rename(ID = id)

dat_4hr <- dat_4hr %>%
  mutate(date = as_datetime(date)) %>%
  rename(ID = id)

dat_8hr <- dat_8hr %>%
  mutate(date = as_datetime(date)) %>%
  rename(ID = id)



# Calculate step lengths and turning angles
dat_1hr_2 <- prepData(dat_1hr, type = 'UTM', coordNames = c('x','y'))
dat_4hr_2 <- prepData(dat_4hr, type = 'UTM', coordNames = c('x','y'))
dat_8hr_2 <- prepData(dat_8hr, type = 'UTM', coordNames = c('x','y'))  #can also be done using lat/long

plot(dat_1hr_2)
plot(dat_4hr_2)
plot(dat_8hr_2)




###############################
### HMMs for 1 hr time step ###
###############################



#### Fit HMM for 2 states using step lengths and turning angles ####

# Calculate displacement separately per ID
dat_1hr_3 <- dat_1hr_2 %>%
  split(.$ID) %>%
  purrr::map(., calc_disp, x, y) %>%
  bind_rows()



## Viz data streams over time per ID

ggplot(dat_1hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_1hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_1hr_3, aes(date, angle)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")



# Pre-define states to set "good" initial values
dat_1hr_3 <- dat_1hr_3 %>%
  mutate(phase = case_when(step > 1 ~ 'Migratory',
                           TRUE ~ 'ARS'))


ggplot(dat_1hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_1hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_1hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_1hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_1hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_1hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(0.2, 3, 0.2, 1) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0.5, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_1hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_1hr_3)
zeromass0 <- c(propzero, 1e-9)  #for zero distances by state
dispPar0 <- c(170, 360, 300, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)


Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_1hr_3 <- dat_1hr_3 %>%
  dplyr::select(-phase)



set.seed(2022)
tic()
fit_1hr_2states <- run.HMMs(data = dat_1hr_3, K = 2, Par0 = Par0,
                            state.names = c('ARS','Migratory'), niter = 20)
toc()  #took 12 min to run


fit_1hr_2states

plot(fit_1hr_2states)
plotPR(fit_1hr_2states, ncores = 5)  #plot of pseudo-residuals look pretty good besides ACF




#### Fit HMM for 3 states using step lengths, turning angles, and displacement ####

# Pre-define states to set "good" initial values
dat_1hr_3 <- dat_1hr_3 %>%
  mutate(phase = case_when(disp < 6 ~ 'Breeding',
                           step > 1 ~ 'Migratory',
                           disp > 6 & step < 1 ~ 'Foraging'))


# Check if a priori assumptions provide "good" rough estimates to set initial values
ggplot(dat_1hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_1hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_1hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_1hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_1hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_1hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(0.15, 0.3, 3, 0.2, 0.2, 1) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0, 0.5, 0.7, 0.99) # (mean_1, mean_2, mean_3, conc_1, conc_2, conc_3)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_1hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_1hr_3)
zeromass0 <- c(propzero, 1e-9, 1e-9)  #for zero distances by state
dispPar0 <- c(1, 570, 400, 1, 200, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)

Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_1hr_3 <- dat_1hr_3 %>%
  dplyr::select(-phase)



set.seed(2022)
tic()
fit_1hr_3states <- run.HMMs(data = dat_1hr_3, K = 3, Par0 = Par0,
                            state.names = c('Breeding','Foraging','Migratory'), niter = 20)
toc()  #took 45 min to run


fit_1hr_3states

plot(fit_1hr_3states)
plotStates(fit_1hr_3states)
timeInStates(fit_1hr_3states)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_1hr_3states, ncores = 5)  #look decent for SL and TA, but Disp could be improved
# results look okay, but not great; qqplot for Disp could be better and ACF plots for SL and Disp show high corr



#### Compare between models ####

AIC(fit_1hr_2states, fit_1hr_3states)  #3 state model is better; consistent w/ state-dependent distributions and mapped states


dat_1hr_3$state <- viterbi(fit_1hr_3states)
dat_1hr_3$state <- factor(dat_1hr_3$state, levels = c(1,3,2))
levels(dat_1hr_3$state) <- c("Breeding","Migratory","Foraging")




###############################
### HMMs for 4 hr time step ###
###############################


#### Fit HMM for 2 states using step lengths, turning angles, and displacement ####

# Calculate displacement separately per ID
dat_4hr_3 <- dat_4hr_2 %>%
  split(.$ID) %>%
  purrr::map(., calc_disp, x, y) %>%
  bind_rows()



## Viz data streams over time per ID

ggplot(dat_4hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_4hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_4hr_3, aes(date, angle)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")



# Pre-define states to set "good" initial values
dat_4hr_3 <- dat_4hr_3 %>%
  mutate(phase = case_when(step > 5 ~ 'Migratory',
                           TRUE ~ 'ARS'))


ggplot(dat_4hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_4hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_4hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_4hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_4hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_4hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(1, 10, 1, 5) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0.4, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_4hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_4hr_3)
zeromass0 <- c(propzero, 1e-9)        #for zero distances by state
dispPar0 <- c(170, 400, 300, 250, zeromass0) # (mu_1, mu_2, sd_1, sd_2, proportion of zeroes likely present per state)

Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_4hr_3 <- dat_4hr_3 %>%
  dplyr::select(-phase)



set.seed(2022)
tic()
fit_4hr_2states <- run.HMMs(data = dat_4hr_3, K = 2, Par0 = Par0,
                            state.names = c('ARS','Migratory'), niter = 20)
toc()  #took 3.5 min to run

fit_4hr_2states

plot(fit_4hr_2states)
plotStates(fit_4hr_2states)
timeInStates(fit_4hr_2states)  #66% ARS, 34% migratory
plotPR(fit_4hr_2states, ncores = 5)  # results look okay, but not great; qqplot for Disp could be better and ACF plots a little better than w/ 1 hr time step tracks






#### Fit HMM for 3 states using step lengths, turning angles, and displacement ####

## Viz data streams over time per ID

ggplot(dat_4hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_4hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_4hr_3, aes(date, angle)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")



# Pre-define states to set "good" initial values
dat_4hr_3 <- dat_4hr_3 %>%
  mutate(phase = case_when(disp < 6 ~ 'Breeding',
                           step > 5 ~ 'Migratory',
                           disp > 6 & step < 5 ~ 'Foraging'))


# Check if a priori assumptions provide "good" rough estimates to set initial values
ggplot(dat_4hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_4hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_4hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_4hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_4hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_4hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(0.5, 1, 12, 0.5, 1, 5) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0, 0.3, 0.6, 0.99) # (mean_1, mean_2, mean_3, conc_1, conc_2, conc_3)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_4hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_4hr_3)
zeromass0 <- c(propzero, 1e-9, 1e-9)  #for zero distances by state
dispPar0 <- c(1, 570, 300, 1, 200, 200, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)

Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_4hr_3 <- dat_4hr_3 %>%
  dplyr::select(-phase)



set.seed(2022)
tic()
fit_4hr_3states <- run.HMMs(data = dat_4hr_3, K = 3, Par0 = Par0,
                            state.names = c('Breeding','Foraging','Migratory'), niter = 20)
toc()  #took 10 min to run

fit_4hr_3states

plot(fit_4hr_3states)
plotStates(fit_4hr_3states)
timeInStates(fit_4hr_3states)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_4hr_3states, ncores = 5)  # results look okay, but not great; qqplot for Disp could be better and ACF plots a little better than w/ 1 hr time step tracks


#### Compare between models ####

AIC(fit_4hr_2states, fit_4hr_3states)  #3 state model is better; consistent w/ state-dependent distributions and mapped states


dat_4hr_3$state <- viterbi(fit_4hr_3states)
dat_4hr_3$state <- factor(dat_4hr_3$state, levels = c(1,3,2))
levels(dat_4hr_3$state) <- c("Breeding","Migratory","Foraging")




###############################
### HMMs for 8 hr time step ###
###############################


#### Fit HMM for 2 states using step lengths, turning angles, and displacement ####

# Calculate displacement separately per ID
dat_8hr_3 <- dat_8hr_2 %>%
  split(.$ID) %>%
  purrr::map(., calc_disp, x, y) %>%
  bind_rows()



## Viz data streams over time per ID

ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_8hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_8hr_3, aes(date, angle)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")



# Pre-define states to set "good" initial values
dat_8hr_3 <- dat_8hr_3 %>%
  mutate(phase = case_when(step > 5 ~ 'Migratory',
                           TRUE ~ 'ARS'))


ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_8hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_8hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_8hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(1, 20, 1, 10) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(3.1, 0, 0.5, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_8hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_8hr_3)
zeromass0 <- c(propzero, 1e-9)        #for zero distances by state
dispPar0 <- c(170, 400, 300, 250, zeromass0) # (mu_1, mu_2, sd_1, sd_2, proportion of zeroes likely present per state)

Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_8hr_3 <- dat_8hr_3 %>%
  dplyr::select(-phase)



set.seed(2022)
tic()
fit_8hr_2states <- run.HMMs(data = dat_8hr_3, K = 2, Par0 = Par0,
                            state.names = c('ARS','Migratory'), niter = 20)
toc()  #took 1.5 min to run

fit_8hr_2states

plot(fit_8hr_2states)
plotStates(fit_8hr_2states)
timeInStates(fit_8hr_2states)  #66% ARS, 34% migratory
plotPR(fit_8hr_2states, ncores = 5)  # results look okay, but not great; qqplot for Disp could be better and ACF plots a little better than w/ 1 hr time step tracks





#### Fit HMM for 3 states using step lengths, turning angles, and displacement ####

# Calculate displacement separately per ID
dat_8hr_3 <- dat_8hr_2 %>%
  split(.$ID) %>%
  purrr::map(., calc_disp, x, y) %>%
  bind_rows()



## Viz data streams over time per ID

ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_8hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_8hr_3, aes(date, angle)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")



# Pre-define states to set "good" initial values
dat_8hr_3 <- dat_8hr_3 %>%
  mutate(phase = case_when(disp < 6 ~ 'Breeding',
                           step > 5 ~ 'Migratory',
                           disp > 6 & step < 5 ~ 'Foraging'))


ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = phase)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_8hr_3, aes(disp, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(step, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(angle, fill = phase)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_8hr_3 %>%
  group_by(phase) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
sum(dat_8hr_3$step == 0)  #check if any SL equal to 0
stepPar0 <- c(1, 1.5, 20, 1, 1, 10) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(3.1, 3.1, 0, 0.4, 0.2, 0.99) # (mean_1, mean_2, mean_3, concentration_1, concentration_2, concentration_3)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_8hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_8hr_3)
zeromass0 <- c(propzero, 1e-9, 1e-9)        #for zero distances by state
dispPar0 <- c(1, 570, 450, 1, 200, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)

Par0 <- list(step = stepPar0, angle = anglePar0, disp = dispPar0)

# Remove 'phase' column to directly compare models by AIC
dat_8hr_3 <- dat_8hr_3 %>%
  dplyr::select(-phase)


set.seed(2022)
tic()
fit_8hr_3states <- run.HMMs(data = dat_8hr_3, K = 3, Par0 = Par0,
                            state.names = c('Breeding','Foraging','Migratory'), niter = 20)
toc()  #took 5 min to run

fit_8hr_3states

plot(fit_8hr_3states)
plotStates(fit_8hr_3states)
timeInStates(fit_8hr_3states)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_8hr_3states, ncores = 5)  # results look okay, but not great; qqplot for Disp could be better and ACF plots a little better than w/ 1 hr time step tracks



#### Compare between models ####

AIC(fit_8hr_2states, fit_8hr_3states)  #3 state model is better; consistent w/ state-dependent distributions and mapped states


# Add state estimate to dataset
dat_8hr_3$state <- viterbi(fit_8hr_3states)
dat_8hr_3$state <- factor(dat_8hr_3$state, levels = c(1,3,2))
levels(dat_8hr_3$state) <- c("Breeding","Migratory","Foraging")






############################
### Fig 3 for manuscript ###
############################

# load spatial layer of brazil
brazil<- ne_countries(scale = 10, country = "Brazil", returnclass = 'sf')

# generate sequence for x axis of density functions
step.seq <- seq(0, max(dat_8hr_3$step, na.rm = TRUE), length = 500)
angle.seq <- seq(-pi, pi, length = 500)
disp.seq <- seq(0, max(dat_8hr_3$disp, na.rm = TRUE), length = 500)

# functions to calculate shape and scale of the gamma distributions from mean and sd
sh <- function(mean, sd) { return(mean^2 / sd^2)}
sc <- function(mean, sd) { return(sd^2 / mean)}

# get density functions of the distributions
step_state.dep.dist <- apply(fit_8hr_3states$mle$step, 2,
                             function(x) dgamma(step.seq, shape = sh(x[1], x[2]),
                                                scale = sc(x[1], x[2]))
) %>%
  cbind(., step.seq) %>%
  data.frame() %>%
  pivot_longer(cols = -step.seq, names_to = "state", values_to = "density") %>%
  rename(x = step.seq) %>%
  mutate(var = 'Step Length (km)')

angle_state.dep.dist <- apply(fit_8hr_3states$mle$angle, 2,
                              function(x) circular::dwrappedcauchy(angle.seq, mu = x[1], rho = x[2])
) %>%
  cbind(., angle.seq) %>%
  data.frame() %>%
  pivot_longer(cols = -angle.seq, names_to = "state", values_to = "density") %>%
  rename(x = angle.seq) %>%
  mutate(var = 'Turning Angle (rad)')

disp_state.dep.dist <- apply(fit_8hr_3states$mle$disp, 2,
                             function(x) dgamma(disp.seq, shape = sh(x[1], x[2]),
                                                scale = sc(x[1], x[2]))
) %>%
  cbind(., disp.seq) %>%
  data.frame() %>%
  pivot_longer(cols = -disp.seq, names_to = "state", values_to = "density") %>%
  rename(x = disp.seq) %>%
  mutate(var = 'Displacement (km)')


state.dep.dist <- rbind(step_state.dep.dist, angle_state.dep.dist, disp_state.dep.dist)
state.dep.dist$var <- factor(state.dep.dist$var, levels = c("Step Length (km)", "Turning Angle (rad)",
                                                            "Displacement (km)"))


# plot state-dependent distributions
state.dep.plot <- ggplot() +
  geom_line(data = state.dep.dist, aes(x = x, y = density, colour = state), size=1) +
  scale_color_manual('', values = MetPalettes$Egypt[[1]][c(1,3,4)]) +
  labs(x = "", y = "Density") +
  theme_bw() +
  theme(legend.text = element_text(size = 14),
        legend.position = "top",
        panel.grid = element_blank(),
        axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 16),
        strip.background = element_blank(),
        strip.placement = "outside") +
  facet_wrap(~ var, scales = "free", strip.position = "bottom")



# plot time-series of state estimates per ID
all.fits <- list(`1 hr` = dat_1hr_3,
                 `4 hr` = dat_4hr_3,
                 `8 hr` = dat_8hr_3) %>%
  bind_rows(.id = 'time.step')


behav.ts.plot <- ggplot() +
  geom_point(data = all.fits %>%
               filter(ID %in% c(205540, 41614)), aes(date, state, color = state), alpha = 0.3) +
  scale_color_manual('', values = MetPalettes$Egypt[[1]][c(1,4,3)], guide = "none") +
  theme_bw() +
  labs(x = "Date", y = "State") +
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  facet_grid(time.step ~ ID, scales = "free_x")



behav.map <- ggplot() +
  geom_sf(data = brazil, fill = "grey60", size = 0.3, color = "black") +
  geom_path(data = all.fits %>%
              filter(ID %in% c(205540, 41614)), aes(lon, lat), alpha = 0.7) +
  geom_point(data = all.fits %>%
               filter(ID %in% c(205540, 41614)), aes(lon, lat, color = state)) +
  scale_color_manual('', values = MetPalettes$Egypt[[1]][c(1,4,3)], guide = "none") +
  coord_sf(xlim = c(-40, -32), ylim = c(-7, -3)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 16)) +
  facet_grid(time.step ~ ID)


state.dep.plot / (behav.ts.plot + behav.map) + plot_annotation(tag_levels = 'a', tag_suffix = ')')

# ggsave("Figures/Fig 3.png", width = 10, height = 7, units = "in", dpi = 400)



#### Export datasets for easy loading ####

# save(dat_1hr_3, dat_4hr_3, dat_8hr_3, fit_1hr_3states, fit_4hr_3states, fit_8hr_3states,
#      file = "Processed_data/HMM_data_and_model_fits.RData")
