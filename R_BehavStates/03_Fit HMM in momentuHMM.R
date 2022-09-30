
#############################################
### Fit discrete-time hidden Markov model ###
#############################################

library(tidyverse)
library(lubridate)
library(momentuHMM)  #v1.5.4
library(sf)  #v1.0.7
library(tictoc)
library(plotly)


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


#### Fit HMM for 3 states using step lengths, turning angles, and displacement ####

# Calculate displacement manually
calc_disp <- function(data, x, y) {
  data$disp <- sqrt((data[,"x"] - data[,"x"][1])^2 + (data[,"y"] - data[,"y"][1])^2)

  return(data)
}

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
  mutate(phase = case_when(disp < 6 ~ 'Breeding',
                           step > 2 ~ 'Migratory',
                           disp > 6 & step < 2 ~ 'Foraging'))


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
stepPar0 <- c(0.2, 0.5, 3.5, 0.2, 0.5, 1) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 0, 0.3, 0.6, 0.99) # (mean_1, mean_2, mean_3, conc_1, conc_2, conc_3)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_1hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_1hr_3)
zeromass0 <- c(propzero, 1e-9, 1e-9)  #for zero distances by state
dispPar0 <- c(1, 570, 350, 1, 200, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)



set.seed(2022)
tic()
fit_1hr_3states <- fitHMM(data = dat_1hr_3,
                          nbStates = 3,
                          dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
                          Par0 = list(step = stepPar0, angle = anglePar0, disp = dispPar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Breeding','Foraging','Migratory'),
                          retryFits = 30)
toc()  #took 4.85 hrs to run

fit_1hr_3states

plot(fit_1hr_3states)
plotStates(fit_1hr_3states)
timeInStates(fit_1hr_3states)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_1hr_3states, ncores = 5)  #look decent for SL and TA, but Disp could be improved
# results look terrible and don't make sense




###############################
### HMMs for 8 hr time step ###
###############################


#### Fit HMM for 3 states using step lengths, turning angles, and displacement ####

# Calculate displacement manually
calc_disp <- function(data, x, y) {
  data$disp <- sqrt((data[,"x"] - data[,"x"][1])^2 + (data[,"y"] - data[,"y"][1])^2)

  return(data)
}

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
                           step > 4 ~ 'Migratory',
                           disp > 6 & step < 4 ~ 'Foraging'))


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
stepPar0 <- c(0.5, 1.5, 15, 1, 1.5, 12) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(3.1, 0, 0, 0.8, 0.1, 0.99) # (mean_1, mean_2, mean_3, concentration_1, concentration_2, concentration_3)

#initial displacement distribution natural scale parameters
whichzero <- which(dat_8hr_3$disp == 0)
propzero <- length(whichzero)/nrow(dat_8hr_3)
zeromass0 <- c(propzero, 1e-9, 1e-9)        #for zero distances by state
dispPar0 <- c(1, 570, 450, 1, 200, 250, zeromass0) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3, proportion of zeroes likely present per state)



set.seed(123)
tic()
fit_8hr_3states <- fitHMM(data = dat_8hr_3,
                                  nbStates = 3,
                                  dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
                                  Par0 = list(step = stepPar0, angle = anglePar0, disp = dispPar0),
                                  formula = ~ 1,
                                  estAngleMean = list(angle=TRUE),
                                  stateNames = c('Breeding','Foraging','Migratory'),
                                  retryFits = 30)
toc()  #took 37.5 min to run

fit_8hr_3states

plot(fit_8hr_3states)
plotStates(fit_8hr_3states)
timeInStates(fit_8hr_3states)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_8hr_3states, ncores = 5)  #look decent for SL and TA, but Disp could be improved




#### Compare between models ####

AIC(fit_hmm_3states_3vars, fit_hmm_3states_3vars2)
AICweights(fit_hmm_3states_3vars, fit_hmm_3states_3vars2)
#AIC actually appears to favor model w/o covariates




#### Export datasets for easy loading ####

save(dat_8hr_2, dat_8hr_3, fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1,
     fit_hmm_3states_covar1, fit_hmm_3states_3vars, fit_hmm_3states_3vars2,
     file = "Processed_data/HMM_data_and_model_fits.RData")
