
#############################################
### Fit discrete-time hidden Markov model ###
#############################################

library(tidyverse)
library(lubridate)
library(momentuHMM)  #v1.5.4
library(sf)  #v1.0.7
library(tictoc)
library(plotly)
library(rerddapXtracto)


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
dat_8hr <- dat_8hr %>%
  mutate(date = as_datetime(date)) %>%
  rename(ID = id)

dat_8hr_2 <- prepData(dat_8hr, type = 'UTM', coordNames = c('x','y'))  #can also be done using lat/long
plot(dat_8hr_2)


# Plot time series and distributions of step lengths
ggplot(dat_8hr_2, aes(step)) +
  geom_histogram(binwidth = 0.5) +
  theme_bw()

ggplot(dat_8hr_2, aes(date, step)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#for 2 states, initial mean values probably around 1 and 30 km whereas SD probably around 0.5 and 5


# Plot time series and distributions of turning angles
ggplot(dat_8hr_2, aes(angle)) +
  geom_histogram(binwidth = pi/8) +
  theme_bw()

ggplot(dat_8hr_2, aes(date, angle)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#for 2 states, initial mean values probably both around 0 rad whereas concentration probably around 0.99 and 0.2





#### Extract sea surface temperature ####

#Available at url: https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.html

# Extract monthly SST
xpos <- dat_8hr_2$lon
ypos <- dat_8hr_2$lat
tpos <- dat_8hr_2$date
sstInfo <- rerddap::info('jplMURSST41')

tic()
sst <- rxtracto(sstInfo, parameter = 'analysed_sst',
                xcoord = xpos, ycoord = ypos, tcoord = tpos,
                xlen = 0.2, ylen = 0.2, progress_bar = TRUE)
toc()  #takes 17 min to run

plotTrack(sst, xpos, ypos, tpos, plotColor = 'thermal')


# Add SST, year-day (yday), and hour to dataset as potential covariates

dat_8hr_2$sst <- sst$`mean analysed_sst`
dat_8hr_2$yday <- yday(dat_8hr_2$date)
dat_8hr_2$hr <- hour(dat_8hr_2$date)


# Is there a relationship between SST and yday?

plot(dat_8hr_2$yday, dat_8hr_2$sst)
#a non-linear pattern does appear to be present; let's stick w/ SST as driver and hr as a cyclical covariate






#### Fit HMM for 2 states using step lengths and turning angles ####

# initial step distribution natural scale parameters
stepPar0 <- c(1, 30, 0.5, 5) # (mu_1, mu_2, sd_1, sd_2)

# initial angle distribution natural scale parameters
anglePar0 <- c(3.1, 0, 0.8, 0.99) # (mean_1, mean_2, concentration_1, concentration_2)


set.seed(123)  #if states get flipped, try adjusting initial params or random seed number (e.g., seed 2022 gave problems w/ flipped states)
tic()
fit_hmm_2states <- fitHMM(data = dat_8hr_2,
                          nbStates = 2,
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Encamped', 'Migratory'),
                          retryFits = 30)  #may be necessary to run more fits
toc()  #took 3.5 min to run

fit_hmm_2states
plot(fit_hmm_2states)
plotPR(fit_hmm_2states, ncores = 5)  #plot of pseudo-residuals show that there are likely some problems; but not terrible





#### Fit HMM for 3 states using step lengths and turning angles ####
# There might be a difference between 'encamped/breeding' behavior and 'foraging' behavior after migration

# Step lengths
ggplot(dat_8hr_2, aes(date, step)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#try means at 1, 3, and 20 km; SDs at 0.5, 2, and 5 km

# Turning angles
ggplot(dat_8hr_2, aes(date, angle)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")
#try means at 0 rad; concentrations at 0.1, 0.5, 0.99



# initial step distribution natural scale parameters
stepPar0 <- c(1, 3, 20, 0.5, 2, 5) # (mu_1, mu_2, mu_3, sd_1, sd_2, sd_3)

# initial angle distribution natural scale parameters
anglePar0 <- c(3.1, 0, 0, 0.8, 0.5, 0.99) # (mean_1, mean_2, mean_3, concentration_1, concentration_2, concentration_3)

set.seed(123)
tic()
fit_hmm_3states <- fitHMM(data = dat_8hr_2,
                          nbStates = 3,
                          dist = list(step = "gamma", angle = "wrpcauchy"),
                          Par0 = list(step = stepPar0, angle = anglePar0),
                          formula = ~ 1,
                          estAngleMean = list(angle=TRUE),
                          stateNames = c('Breeding', 'Foraging', 'Migratory'),
                          retryFits = 30)
toc()  #took 16.5 min to run

fit_hmm_3states

plot(fit_hmm_3states)  #model can't really differentiate between the breeding and foraging ground movements; but does appear to detect a 3rd state
plotStates(fit_hmm_3states)  #possible that 1st state is 'Encamped' and 2nd state is 'ARS'
timeInStates(fit_hmm_3states)  #primarily in Encamped/Breeding or Foraging state; only 3% Migratory
plotPR(fit_hmm_3states, ncores = 5)  #ACF and QQ plot for SL looks much better





#### Fit HMM for 2 states including covariates ####


# formula for transition probabilities
formula <- ~ sst * cosinor(hr, period = 24)

# initial parameters (obtained from nested model first model)
Par0_covar <- getPar0(model = fit_hmm_2states, formula = formula)


set.seed(123)
tic()
fit_hmm_2states_covar1 <- fitHMM(data = dat_8hr_2,
                                nbStates = 2,
                                dist = list(step = "gamma", angle = "wrpcauchy"),
                                Par0 = Par0_covar$Par,
                                beta0 = Par0_covar$beta,
                                formula = formula,
                                estAngleMean = list(angle=TRUE),
                                stateNames = c('Encamped', 'Migratory'),
                                retryFits = 30)
toc()  #took 6 min to run

fit_hmm_2states_covar1
plot(fit_hmm_2states_covar1, plotCI = TRUE, plotStationary = TRUE)  #no noticeable impact of including SST, although slight effect of yday; stationary states show some patterns
plotStates(fit_hmm_2states_covar1)
stationary(fit_hmm_2states_covar1)  #gets stationary state probs plotted above per obs
timeInStates(fit_hmm_2states_covar1)  #95% encamped; 5% migratory
plotPR(fit_hmm_2states_covar1, ncores = 5)  #looks about the same as earlier 2-state model









#### Fit HMM for 3 states including covariates ####


# formula for transition probabilities
formula <- ~ sst * cosinor(hr, period = 24)

# initial parameters (obtained from nested model first model)
Par0_covar3 <- getPar0(model = fit_hmm_3states, formula = formula)


set.seed(123)
tic()
fit_hmm_3states_covar1 <- fitHMM(data = dat_8hr_2,
                                 nbStates = 3,
                                 dist = list(step = "gamma", angle = "wrpcauchy"),
                                 Par0 = Par0_covar3$Par,
                                 beta0 = Par0_covar3$beta,
                                 formula = formula,
                                 estAngleMean = list(angle=TRUE),
                                 stateNames = c('Breeding', 'Foraging', 'Migratory'),
                                 retryFits = 30)
toc()  #took 25.5 min to run

fit_hmm_3states_covar1
plot(fit_hmm_3states_covar1, plotCI = TRUE, plotStationary = TRUE)  #seems to also perform worse
plotStates(fit_hmm_3states_covar1)
stationary(fit_hmm_3states_covar1)  #gets stationary state probs plotted above per obs
timeInStates(fit_hmm_3states_covar1)  #primarily Breeding or Foraging; 5% migratory
plotPR(fit_hmm_3states_covar1, ncores = 5)  #looks like one of best model diagnostics so far





### Compare among models ###

AIC(fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1, fit_hmm_3states_covar1)
AICweights(fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1, fit_hmm_3states_covar1)
# of these models, AIC suggests that the most complicated one is best; but not necessarily the case if not interpretable





#### Fit HMM for 3 states that evaluates SL, TA, and NSD ####

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



## Viz displacement and step length over time per ID

ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")

ggplot(dat_8hr_3, aes(date, step)) +
  geom_path(aes(group = ID)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free")


# Pre-define states to set "good" initial values
dat_8hr_3 <- dat_8hr_3 %>%
  mutate(state = case_when(disp < 6 ~ 'Breeding',
                           step > 4 ~ 'Migratory',
                           disp > 6 & step < 4 ~ 'Foraging'))


ggplot(dat_8hr_3, aes(date, disp)) +
  geom_path(aes(group = ID, color = state)) +
  theme_bw() +
  facet_wrap(~ID, scales = "free_x")


ggplot(dat_8hr_3, aes(disp, fill = state)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(step, fill = state)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()

ggplot(dat_8hr_3, aes(angle, fill = state)) +
  geom_histogram(alpha = 0.6) +
  theme_bw()


dat_8hr_3 %>%
  group_by(state) %>%
  summarize(mean.step = mean(step, na.rm = T),
            sd.step = sd(step, na.rm = T),
            mean.disp = mean(disp, na.rm = T),
            sd.disp = sd(disp, na.rm = T))




# initial step distribution natural scale parameters
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
fit_hmm_3states_3vars <- fitHMM(data = dat_8hr_3,
                                  nbStates = 3,
                                  dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
                                  Par0 = list(step = stepPar0, angle = anglePar0, disp = dispPar0),
                                  formula = ~ 1,
                                  estAngleMean = list(angle=TRUE),
                                  stateNames = c('Breeding','Foraging','Migratory'),
                                  retryFits = 30)
toc()  #took 37.5 min to run

fit_hmm_3states_3vars

plot(fit_hmm_3states_3vars)
plotStates(fit_hmm_3states_3vars)
timeInStates(fit_hmm_3states_3vars)  #66% breeding, 29% foraging, 5% migratory
plotPR(fit_hmm_3states_3vars, ncores = 5)  #look decent for SL and TA, but Disp could be improved







#### Fit HMM for 3 states that evaluates SL, TA, and NSD w/ covariates ####


# initial parameters (obtained from nested model; i.e., previous model)
Par0_covar4 <- getPar0(model = fit_hmm_3states_3vars, formula = ~ sst)


set.seed(123)
tic()
fit_hmm_3states_3vars2 <- fitHMM(data = dat_8hr_3,
                                 nbStates = 3,
                                 dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
                                 Par0 = Par0_covar4$Par,
                                 beta0 = Par0_covar4$beta,
                                 formula = ~ sst,
                                 # DM = DM,
                                 estAngleMean = list(angle=TRUE),
                                 stateNames = c('Breeding','Foraging','Migratory'),
                                 retryFits = 30)
toc()  #took 32 min to run

fit_hmm_3states_3vars2
plot(fit_hmm_3states_3vars2, plotCI = TRUE, plotStationary = TRUE)  #affects of covariates are greater on SL and TA vs TPM; stationary states show some patterns
plotStates(fit_hmm_3states_3vars2)
stationary(fit_hmm_3states_3vars2)  #gets stationary state probs plotted above per obs
timeInStates(fit_hmm_3states_3vars2)
plotPR(fit_hmm_3states_3vars2, ncores = 5)  #looks about the same as before





#### Compare between models ####

AIC(fit_hmm_3states_3vars, fit_hmm_3states_3vars2)
AICweights(fit_hmm_3states_3vars, fit_hmm_3states_3vars2)
#AIC actually appears to favor model w/o covariates




#### Export datasets for easy loading ####

save(dat_8hr_2, dat_8hr_3, fit_hmm_2states, fit_hmm_3states, fit_hmm_2states_covar1,
     fit_hmm_3states_covar1, fit_hmm_3states_3vars, fit_hmm_3states_3vars2,
     file = "Processed_data/HMM_data_and_model_fits.RData")
