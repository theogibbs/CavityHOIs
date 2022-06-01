## Packages

library(deSolve)
library(MASS)
library(tidyverse)
library(plyr)
library(reshape2)
library(nleqslv)
library(VGAM)
library(truncnorm)
library(gridExtra)

## Functions

# multiplies matrices so we can use apply easily later
MatMult <- function(M, x) return(M %*% x)

# removes the diagonal elements of the B tensor that correspond to cubic self-regulation
NoSelfHOIs <- function(M) {
  newM <- M
  for(i in 1:dim(newM)[3]) {
    newM[i,i,i] <- 0
  }
  return(newM)
}

# creates the pairwise interaction matrix A with desired statistics
BuildA <- function(S = 50, mu_A = -2, sigma_A = 0.5, rho_A = 0) {
  
  NumPairs <- S * (S - 1) / 2
  scaled_mu_A <- mu_A / S
  mus <- c(scaled_mu_A, scaled_mu_A)
  scaled_sigma_A <- sigma_A / sqrt(S)
  
  covariance.matrix <- matrix(c(scaled_sigma_A^2,
                                rho_A * scaled_sigma_A^2,
                                rho_A * scaled_sigma_A^2,
                                scaled_sigma_A^2),
                              2, 2)
  Pairs <- mvrnorm(NumPairs, mus, covariance.matrix)
  
  A <- matrix(0, S, S)
  k <- 1
  for (i in 1:(S-1)){
    for (j in (i + 1):S){
      A[j,i] <- Pairs[k,1]
      A[i,j] <- Pairs[k,2]
      k <- k + 1
    }
  }
  return(A)
}

# creates the HOI tensor with desired statistics and correlation to the pairwise interaction matrix A
BuildB <- function(mu_B = -2, sigma_B = 0.5, rho_B = 0, mu_A = -2, A) {
  
  S <- nrow(A)
  scaled_mu_B <- mu_B / S^2
  scaled_sigma_B <- sigma_B / S

  B <- array(rnorm(S*S*S, mean = scaled_mu_B, sd = scaled_sigma_B), c(S, S, S))
  centered_A <- A - matrix(mu_A / S, nrow = S, ncol = S)
  diag(centered_A) <- 0
  B <- B + rho_B * array(centered_A, c(S, S, S)) / S
  
  B <- NoSelfHOIs(B)
  return(B)
}

# creates a list of parameters for integrating the dynamics
BuildPars <- function(input_params) {
  
  set.seed(input_params$ParsID)
  pars <- with(input_params, list(S = S, r = rnorm(S, mean = MuR, sd = SigmaR),
                                  d = rnorm(S, mean = MuD, sd = SigmaD),
                                  A = BuildA(S, MuA, SigmaA, RhoA),
                                  h = h))
  pars$B <- with(input_params, BuildB(MuB, SigmaB, RhoB, MuA, pars$A))

  return(pars)
}

# SUPPLEMENTARY: creates a list of parameters for integrating the dynamics witout non-normal distributions
BuildNonNormalPars <- function(input_params) {
  
  set.seed(input_params$ParsID)
  pars <- with(input_params, list(S = S, r = runif(S, min = MuR - sqrt(3) * SigmaR, max = MuR + sqrt(3) * SigmaR),
                                  d = runif(S, min = MuD - sqrt(3) * SigmaD, max = MuD + sqrt(3) * SigmaD),
                                  h = h))
  pars$A <- with(input_params, {
    scaled_mu_A <- MuA / S
    scaled_sigma_A <- SigmaA / sqrt(S)
    A <- matrix(runif(S^2, min = scaled_mu_A - sqrt(3) * scaled_sigma_A, max = scaled_mu_A + sqrt(3) * scaled_sigma_A), S, S)
    return(A)
    })
  pars$B <- with(input_params, {
    scaled_mu_B <- MuB / S^2
    scaled_sigma_B <- SigmaB / S
    
    B <- array(runif(S*S*S, min = scaled_mu_B - sqrt(3) * scaled_sigma_B, max = scaled_mu_B + sqrt(3) * scaled_sigma_B), c(S, S, S))
    B <- NoSelfHOIs(B)
    return(B)
    })
  
  return(pars)
}

# computes the growth rates of the current state and parameters for the linear ODEs
GetGrowthRates <- function(N, pars) {
  
  S <- pars$S
  r <- pars$r
  d <- pars$d
  A <- pars$A
  B <- pars$B
  h <- pars$h
  
  B <- matrix(as.vector(B), nrow = S, ncol = S^2)
  Ncombs <- as.vector(outer(N, N))
  growth_rates <- r - d * N + A %*% N + B %*% (Ncombs / (1 + h * Ncombs))
  
  return(growth_rates)
}

# SUPPLEMENTARY: computes the growth rates of the current state and parameters with cubic self-regulation
GetCubicGrowthRates <- function(N, pars) {
  
  S <- pars$S
  r <- pars$r
  d <- pars$d
  A <- pars$A
  B <- pars$B
  h <- pars$h
  
  B <- matrix(as.vector(B), nrow = S, ncol = S^2)
  Ncombs <- as.vector(outer(N, N))
  growth_rates <- r - d * N^2 + A %*% N + B %*% (Ncombs / (1 + h * Ncombs))
  
  return(growth_rates)
}

# SUPPLEMENTARY: computes the growth rates of the current state and parameters for the min-max dynamics
GetMinMaxGrowthRates <- function(N, pars) {
  
  S <- pars$S
  r <- pars$r
  d <- pars$d
  A <- pars$A
  B <- pars$B
  h <- pars$h
  
  B <- matrix(as.vector(B), nrow = S, ncol = S^2)
  Nks <- rep(N, each = S)
  sum_ints <- B %*% matrix(Nks, nrow = S^2, ncol = S)
  Bmat <- sum_ints / (h + sum_ints)
  Bmat[sum_ints < 0] <- 0
  Bmat <- 1 - Bmat
  growth_rates <- r - d * N + (A * Bmat) %*% N
  
  return(growth_rates)
}

# SUPPLEMENTARY: computes the growth rates of the current state and parameters for the min-max dynamics
GetConstrainedGrowthRates <- function(N, pars) {
  
  S <- pars$S
  r <- pars$r
  d <- pars$d
  A <- pars$A
  B <- pars$B
  h <- pars$h
  
  B <- matrix(as.vector(B), nrow = S, ncol = S^2)
  Nks <- rep(N, each = S)
  ho_ints <- B %*% matrix(Nks, nrow = S^2, ncol = S)
  ho_ints <- ho_ints %*% diag(N, S, S)
  
  pw_ints <- A %*% diag(N, S, S)
  
  ints <- pw_ints + ho_ints
  
  pos_ints <- ints[pw_ints > 0]
  pw_pos_ints <- pw_ints[pw_ints > 0]
  pos_ints <- pmax(0, pos_ints)
  pos_ints <- pmin(2 * pw_pos_ints, pos_ints)
  ints[pw_ints > 0] <- pos_ints
  
  neg_ints <- ints[pw_ints < 0]
  pw_neg_ints <- pw_ints[pw_ints < 0]
  neg_ints <- pmin(0, neg_ints)
  neg_ints <- pmax(2 * pw_neg_ints, neg_ints)
  ints[pw_ints < 0] <- neg_ints
  diag(ints) <- 0
  
  growth_rates <- r - d * N + rowSums(ints)
  
  return(growth_rates)
}

# wrapper to compute ODE derivatives
Dynamics <- function(time, state, params) {
  dNdt <- state * GetGrowthRates(N = state, pars = params)
  return(list(dNdt))
}

# SUPPLEMENTARY: wrapper to compute derivatives for cubic self-regulation
CubicDynamics <- function(time, state, params) {
  dNdt <- state * GetCubicGrowthRates(N = state, pars = params)
  return(list(dNdt))
}

# SUPPLEMENTARY: wrapper to compute derivatives for cubic self-regulation
SatCorrDynamics <- function(time, state, params) {
  dNdt <- state * GetMinMaxGrowthRates(N = state, pars = params)
  return(list(dNdt))
}

# wrapper to compute derivatives for cubic self-regulation
ConstrainedDynamics <- function(time, state, params) {
  dNdt <- state * GetConstrainedGrowthRates(N = state, pars = params)
  return(list(dNdt))
}

# integrates the dynamics using deSolve
IntegrateDynamics <- function(inistate, pars, endtime, timestep, fn){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars))  
  return(timeseries)
}

# plots time series
PlotSeries <- function(series, title = "Community Dynamics") {
  meltseries <- melt(series, id.vars = "time")
  pl <- ggplot() +
    geom_line(data = meltseries, aes(x = time, y = value, color = variable),
              size = 3, alpha = 0.5) +
    ggtitle(title) +
    xlab("Time") + ylab("Abundance") +
    theme_classic() +
    theme(legend.position = "none",
          text = element_text(size=30),
          strip.text.x = element_text(size = 25),
          strip.text.y = element_text(size = 25),
          legend.text=element_text(size = 25))
  return(pl)
}

# runs the dynamics, gets the endpoints and checks if the system is at an uninvadeable equilibrium
GetAbds <- function(pars, settings) {
  
  S <- pars$S
  inistate <- with(settings, runif(S, min = inimin, max = inimax))
  out <- with(settings, IntegrateDynamics(inistate, pars, endtime, endtime, Dynamics))
  
  end_state <- as.numeric(out[(nrow(out)),-1])
  survivors <- which(end_state > settings$abd_cutoff)
  invaders <- setdiff(1:S, survivors)
  end_state[invaders] <- 0
  
  growth_rates <- GetGrowthRates(end_state, pars)
  survivor_grs <- growth_rates[survivors]
  equilibrium <- as.logical(prod(survivor_grs <= settings$gr_cutoff))
  
  invader_grs <- growth_rates[invaders]
  pos_invaders <- which(invader_grs > 0)
  uninvadeable <- is_empty(pos_invaders)
  
  ret <- data.frame(SpeciesID = 1:S,
                    Abundances = end_state,
                    Equilibrium = equilibrium,
                    Uninvadeable = uninvadeable)
  return(ret)
}

# SUPPLEMENTARY: same as above function but with cubic self-regulation
GetCubicAbds <- function(pars, settings) {
  
  S <- pars$S
  inistate <- with(settings, runif(S, min = inimin, max = inimax))
  out <- with(settings, IntegrateDynamics(inistate, pars, endtime, endtime, CubicDynamics))
  
  end_state <- as.numeric(out[(nrow(out)),-1])
  survivors <- which(end_state > settings$abd_cutoff)
  invaders <- setdiff(1:S, survivors)
  end_state[invaders] <- 0
  
  growth_rates <- GetCubicGrowthRates(end_state, pars)
  survivor_grs <- growth_rates[survivors]
  equilibrium <- as.logical(prod(survivor_grs <= settings$gr_cutoff))
  
  invader_grs <- growth_rates[invaders]
  pos_invaders <- which(invader_grs > 0)
  uninvadeable <- is_empty(pos_invaders)
  
  ret <- data.frame(SpeciesID = 1:S,
                    Abundances = end_state,
                    Equilibrium = equilibrium,
                    Uninvadeable = uninvadeable)
  return(ret)
}

# SUPPLEMENTARY: same as above function but with correlated and saturating HOIs
GetSatCorrAbds <- function(pars, settings) {
  
  S <- pars$S
  inistate <- with(settings, runif(S, min = inimin, max = inimax))
  out <- with(settings, IntegrateDynamics(inistate, pars, endtime, endtime, SatCorrDynamics))
  
  end_state <- as.numeric(out[(nrow(out)),-1])
  survivors <- which(end_state > settings$abd_cutoff)
  invaders <- setdiff(1:S, survivors)
  end_state[invaders] <- 0
  
  growth_rates <- GetMinMaxGrowthRates(end_state, pars)
  survivor_grs <- growth_rates[survivors]
  equilibrium <- as.logical(prod(survivor_grs <= settings$gr_cutoff))
  
  invader_grs <- growth_rates[invaders]
  pos_invaders <- which(invader_grs > 0)
  uninvadeable <- is_empty(pos_invaders)
  
  ret <- data.frame(SpeciesID = 1:S,
                    Abundances = end_state,
                    Equilibrium = equilibrium,
                    Uninvadeable = uninvadeable)
  return(ret)
}

# SUPPLEMENTARY: same as above function but with HOIs constrained to be less than pairwise interactions
GetConstrainedAbds <- function(pars, settings) {
  
  S <- pars$S
  inistate <- with(settings, runif(S, min = inimin, max = inimax))
  out <- with(settings, IntegrateDynamics(inistate, pars, endtime, endtime, ConstrainedDynamics))
  
  end_state <- as.numeric(out[(nrow(out)),-1])
  survivors <- which(end_state > settings$abd_cutoff)
  invaders <- setdiff(1:S, survivors)
  end_state[invaders] <- 0
  
  growth_rates <- GetConstrainedGrowthRates(end_state, pars)
  survivor_grs <- growth_rates[survivors]
  equilibrium <- as.logical(prod(survivor_grs <= settings$gr_cutoff))
  
  invader_grs <- growth_rates[invaders]
  pos_invaders <- which(invader_grs > 0)
  uninvadeable <- is_empty(pos_invaders)
  
  ret <- data.frame(SpeciesID = 1:S,
                    Abundances = end_state,
                    Equilibrium = equilibrium,
                    Uninvadeable = uninvadeable)
  return(ret)
}

# creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  pars_list <- plyr::alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetAbds, settings)
  
  return(ret_abds)
}

# SUPPLEMENTARY: creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverCubicParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  pars_list <- plyr::alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetCubicAbds, settings)
  
  return(ret_abds)
}

# SUPPLEMENTARY: creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverSatCorrParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  pars_list <- plyr::alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetSatCorrAbds, settings)
  
  return(ret_abds)
}

# SUPPLEMENTARY: creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverNonNormalParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  pars_list <- plyr::alply(.data = iterated_params, .margins = 1, .fun = BuildNonNormalPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetAbds, settings)
  
  return(ret_abds)
}

# SUPPLEMENTARY: creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverConstrainedParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  pars_list <- plyr::alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetConstrainedAbds, settings)
  
  return(ret_abds)
}

# computes the cavity method equations
CavitySoln <- function(x, cur_data) {
  
  phi <- x[1]
  avgN <- x[2]
  secN <- x[3]
  v <- x[4]
  
  avg_norm <- with(cur_data, v / (MuD + MuA / S) * (MuR + phi * MuA * avgN + S * (S - 1) / S^2 * MuB * phi^2 * avgN^2))
  var_norm <- with(cur_data,
                   v^2 / (MuD + MuA / S)^2 * (SigmaR^2 + (S - 1) / S * SigmaA^2 * phi * secN + 2 * phi^2 * avgN * secN * SigmaA^2 * RhoB
                          + S * (S - 1) / S^2 * (SigmaB^2 + RhoB^2 * SigmaA^2) * phi^2 * secN^2))
  
  
  if(var_norm < 1e-6) var_norm <- 1e-6
  error_fn <- erf(avg_norm / sqrt(2 * var_norm))
  
  eq1 <- 0.5 * (1 + error_fn)
  eq2 <- avg_norm / 2 * (1 + error_fn) + sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  eq3 <- (avg_norm^2 + var_norm) / 2 * (1 + error_fn)
  eq3 <- eq3 + avg_norm * sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  
  eq1 <- phi - eq1
  eq2 <- avgN - (1 / phi) * eq2
  eq3 <- secN - (1 / phi) * eq3
  eq4 <- with(cur_data, 1 - v * (MuD - phi * RhoA * SigmaA^2 * v))
  
  return(c(eq1, eq2, eq3, eq4))
}

# solves the cavity method equations
GetPredictions <- function(out_data, max_trials) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredSec = c(), PredV = c(), FuncVal = c(), TermCode = c())
  
  for(row in 1:nrow(out_data)) {
    
    if(row %% 10 == 0) {
      print("Row number:")
      print(row)
    }
    cur_data <- out_data[row,]
    
    ini_phi <- 1
    
    f <- function(ini_mean, cur_data) return(with(cur_data, MuR + (MuA - 1) * ini_mean + MuB * ini_mean^2))
    testxs <- seq(0, 2, length.out = 1000); fxs <- c(); gxs <- c(); hxs <- c()
    for(curx in testxs) fxs <- c(fxs, f(curx, cur_data))
    max_mean <- testxs[fxs < 0][1]
    ini_mean <- uniroot(f, c(0, max_mean), cur_data)$root
    
    g <- function(ini_sec, cur_data) {
      return(with(cur_data, SigmaR^2 + ini_mean^2 + (SigmaA^2 - 1) * ini_sec + SigmaB^2 * ini_sec^2))}
    for(curx in testxs) gxs <- c(gxs, g(curx, cur_data))
    max_sec <- testxs[gxs < 0][1]
    if(!is.na(max_sec)) {
      ini_sec <- uniroot(g, c(0, max_sec), cur_data)$root
    } else {
      ini_sec <- ini_mean
    }
    
    h <- function(ini_v, cur_data) with(cur_data, 1 - ini_v * MuD)
    for(curx in testxs) hxs <- c(hxs, h(curx, cur_data))
    max_v <- testxs[hxs < 0][1]
    ini_v <- uniroot(h, c(0, max_v), cur_data)$root
    
    ini_guess <- c(1 - 0.25 * with(cur_data, sqrt(SigmaA^2 + SigmaB^2)), ini_mean, ini_sec, ini_v)
    cav_soln <- nleqslv(ini_guess, CavitySoln, method = "Broyden", jac = NULL, cur_data)
    fval <- sum(cav_soln$fvec)
    term_code <- cav_soln$termcd
    cav_soln <- cav_soln$x
    
    trial <- 1
    
    while(term_code != 1 && trial < max_trials) {
      ini_guess <- runif(4, 0, 1)
      cav_soln <- nleqslv(ini_guess, CavitySoln, method = "Broyden", jac = NULL, cur_data)
      fval <- sum(cav_soln$fvec)
      term_code <- cav_soln$termcd
      cav_soln <- cav_soln$x
      trial <- trial + 1
    }
    
    PRED <- rbind(PRED, data.frame(PredFraction = ifelse(term_code == 1, cav_soln[1], NA), PredMean = cav_soln[2],
                                   PredSec = cav_soln[3], PredV = cav_soln[4], FuncVal = fval, TermCode = term_code))
  }
  
  ret_out_data <- dplyr::bind_cols(out_data, PRED)
  return(ret_out_data)
}

# checks whether or not abundances are at equilibrium
CheckAbds <- function(abds) {
  ret <- 0
  
  EqRuns <- abds %>%
    filter(Equilibrium == FALSE) %>%
    summarise(NonEqNumRuns = length(Equilibrium))
  
  if(EqRuns$NonEqNumRuns != 0) {
    print("WARNING: Some runs did not reach equilibrium.")
    ret <- 1
  }
  
  InvRuns <- abds %>%
    filter(Equilibrium == TRUE, Uninvadeable == FALSE) %>%
    summarise(InvNumRuns = length(Equilibrium))
  
  if(InvRuns$InvNumRuns != 0) {
    print("WARNING: Some equilibria are invasible.")
    ret <- 1
  }
  return(ret)
}

# labels abundance data by the type of interaction and summarizes their statistics
LabelAbds <- function(abds) {
  
  ret_abds <- abds %>%
    mutate(Interaction= ifelse(SigmaA != 0, ifelse(SigmaB == 0, "Pairwise", "Mixed"), "Higher Order")) %>%
    mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Mixed", "Higher Order"))) %>%
    mutate(Mu = MuA + MuB, Sigma = signif(sqrt(SigmaA^2 + SigmaB^2 + RhoB^2 * SigmaA^2 + RhoB * SigmaA^2)))
  
  return(ret_abds)
}

# computes statistics over replicates
GetStatistics <- function(abds) {
  ret_stats <- abds %>% dplyr::group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, h, CommunityID) %>%
    dplyr::summarise(CommPhi = unique(sum(Abundances > 0) / length(Abundances)), CommMean = mean(Abundances[Abundances > 0]),
                     CommSecAbd = mean(Abundances[Abundances > 0]^2), CommVar = CommSecAbd - CommMean^2,
                     CommFourthAbd = mean(Abundances[Abundances > 0]^4), Interaction = unique(Interaction)) %>%
    dplyr::group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, h) %>%
    dplyr::summarise(Phi = mean(CommPhi), MeanAbd = mean(CommMean), SecAbd = mean(CommSecAbd), VarAbd = mean(CommVar),
                     FourthAbd = mean(CommFourthAbd), ErrorPhi = sd(CommPhi), ErrorMean = sd(CommMean), ErrorSec = sd(CommSecAbd),
                     ErrorVar = sd(CommVar), ErrorFourth = sd(CommFourthAbd), Mu = unique(MuA + MuB),
                     Sigma = unique(signif(sqrt(SigmaA^2 + SigmaB^2))),
                     Interaction = unique(Interaction), NumRepl = length(CommPhi))
  return(ret_stats)
}

# computes the statitsics of the non-truncated normal dist
# from the statistics of the truncated normal
ParentNorm <- function(x, pred_mean, pred_sec) {
  mu <- x[1]
  sigma <- x[2]
  pred_var <- pred_sec - pred_mean^2
  
  mu_eq <- pred_mean - etruncnorm(a = 0, b = Inf, mean = mu, sd = sigma)
  var_eq <- pred_var - vtruncnorm(a = 0, b = Inf, mean = mu, sd = sigma)
  
  return(c(mu_eq, var_eq))
}

# adjusts the mean and variance of the coexisting species to be the mean and variance
# of the non-truncationed normal distribution for the cavity method predictions
GetAdjPreds <- function(pred_data) {
  pred_abds <- pred_data
  pred_abds$AdjPredMean <- 0
  pred_abds$AdjPredSD <- 0
  pred_abds$PredSD <- sqrt(abs(pred_abds$PredSec - pred_abds$PredMean^2))
  
  for(i in 1:nrow(pred_abds)) {
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    new_preds <- nleqslv(c(pred_mean, sqrt(abs(pred_sec - pred_mean^2))), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$AdjPredSD[i] <- new_preds[2]
  }
  
  pred_abds$ExpMin <- pred_abds$PredMean - pred_abds$PredSD * sqrt(2 * log(pred_abds$S))
  pred_abds$AdjExpMin <- pred_abds$AdjPredMean - pred_abds$AdjPredSD * sqrt(2 * log(pred_abds$S))
  
  return(pred_abds)
}

# adjusts the mean and variance of the coexisting species to be the mean and variance
# of the non-truncationed normal distribution for both the cavity method predictions
# and the simulation data
GetAdjStats <- function(pred_data) {
  pred_abds <- pred_data
  pred_abds$AdjMean <- 0
  pred_abds$AdjSD <- 0
  pred_abds$AdjPredMean <- 0
  pred_abds$AdjPredSD <- 0
  
  pred_abds$SD <- sqrt(abs(pred_abds$VarAbd))
  pred_abds$PredSD <- sqrt(abs(pred_abds$PredSec - pred_abds$PredMean^2))
  
  for(i in 1:nrow(pred_abds)) {
    cur_mean <- pred_abds$MeanAbd[i]
    cur_sec <- pred_abds$SecAbd[i]
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    
    new_abds <- nleqslv(c(cur_mean, sqrt(abs(pred_sec - pred_mean^2))), fn = ParentNorm, jac = NULL, cur_mean, cur_sec)$x
    pred_abds$AdjMean[i] <- new_abds[1]
    pred_abds$AdjSD[i] <- new_abds[2]
    
    new_preds <- nleqslv(c(pred_mean, sqrt(abs(pred_sec - pred_mean^2))), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$AdjPredSD[i] <- new_preds[2]
  }
  
  pred_abds$ExpMin <- pred_abds$PredMean - pred_abds$PredSD * sqrt(2 * log(pred_abds$S))
  pred_abds$AdjExpMin <- pred_abds$AdjPredMean - pred_abds$AdjPredSD * sqrt(2 * log(pred_abds$S))
  
  return(pred_abds)
}


