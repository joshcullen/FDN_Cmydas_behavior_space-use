
# Calculate displacement manually
calc_disp <- function(data, x, y) {
  data$disp <- sqrt((data[,"x"] - data[,"x"][1])^2 + (data[,"y"] - data[,"y"][1])^2)

  return(data)
}

#-----------------------------------------


### Function to test different sets of initial values in iterative manner
run.HMMs.internal = function(data, K, Par0, state.names, p, seeds) {

  set.seed(seeds)

      # Step length
      stepMean0 <- runif(K,
                         min = Par0$step[1:K] / 2,
                         max = Par0$step[1:K] * 2)
      stepSD0 <- runif(K,
                       min = Par0$step[(K+1):(K*2)] / 2,
                       max = Par0$step[1:K] * 2)
      whichzero_sl <- which(data$step == 0)
      propzero_sl <- length(whichzero_sl)/nrow(data)
      zeromass0_sl <- c(propzero_sl, rep(0, K-1))        #for zero distances by state


      # Displacement
      dispMean0 <- runif(K,
                         min = Par0$disp[1:K] / 2,
                         max = Par0$disp[1:K] * 2)
      dispSD0 <- runif(K,
                       min = Par0$disp[(K+1):(K*2)] / 2,
                       max = Par0$disp[1:K] * 2)
      whichzero_disp <- which(data$disp == 0)
      propzero_disp <- length(whichzero_disp)/nrow(data)
      zeromass0_disp <- c(propzero_disp, rep(0, K-1))        #for zero distances by state


      # Fit model
      if(propzero_sl > 0) {  #don't include zero mass if no 0s present
        stepPar0 <- c(stepMean0, stepSD0, zeromass0_sl)
      } else {
        stepPar0 <- c(stepMean0, stepSD0)
      }

      anglePar0 <- Par0$angle

      if(propzero_disp > 0) {  #don't include zero mass if no 0s present
        dispPar0 <- c(dispMean0, dispSD0, zeromass0_disp)
      } else {
        dispPar0 <- c(dispMean0, dispSD0)
      }

      hmm.res <- fitHMM(data = data, nbStates = K,
                          Par0 = list(step = stepPar0, angle = anglePar0, disp = dispPar0),
                          dist = list(step = "gamma", angle = "wrpcauchy", disp = "gamma"),
                          formula = ~ 1, stationary=TRUE, #stationary for a slightly better fit
                          estAngleMean = list(angle=TRUE),
                          stateNames = state.names#,
                          # optMethod = "Nelder-Mead"
                          )

      # Update progress bar
      p()

      return(hmm.res)
}


#-----------------------------------


run.HMMs = function(data, K, Par0, state.names, niter) {

  # Convert data.frame into list of identical elements to fit w/ different initial values
  hmm.list <- lapply(1:niter, function(x) data)

  # Fit model
  plan(multisession, workers = availableCores() - 2)
  seeds <- future.apply::future_lapply(seq_along(hmm.list), FUN = function(x) sample(1:1e6, 1),
                                       future.chunk.size = Inf, future.seed = TRUE)  #set seed per list element

  handlers(handler_progress(incomplete=".", complete="*", current="o", clear = FALSE))
  progressr::with_progress({
    #set up progress bar
    p<- progressr::progressor(steps = length(hmm.list))

  hmm.res <- furrr::future_map2(hmm.list, seeds,
                               ~run.HMMs.internal(data = .x,
                                                  K = K,
                                                  Par0 = Par0,
                                                  state.names = state.names,
                                                  p = p,
                                                  seeds = .y),
                               .options = furrr_options(seed = TRUE))
  })

  plan(sequential)

    # Extract likelihoods of fitted models
    allnllk <- unlist(lapply(hmm.res, function(m) m$mod$minimum))

    # Index of best fitting model (smallest negative log-likelihood)
    whichbest <- which.min(allnllk)

    # Best fitting model
    res <- hmm.res[[whichbest]]

  return(res)
}

