## Packages

library(deSolve)
library(MASS)
library(plyr)
library(tidyverse)
library(reshape2)

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

  pars <- with(input_params, list(S = S, r = rnorm(S, mean = MuR, sd = SigmaR),
                                  d = rnorm(S, mean = MuD, sd = SigmaD),
                                  A = BuildA(S, MuA, SigmaA, RhoA)))
  pars$B <- with(input_params, BuildB(MuB, SigmaB, RhoB, MuA, pars$A))
  
  return(pars)
}


# computes the growth rates of the current state and parameters
GetGrowthRates <- function(N, pars) {
  
  r <- pars$r
  d <- pars$d
  A <- pars$A
  B <- pars$B
  
  M <- A + apply(B, 3, MatMult, N)
  
  growth_rates <- (r - d * N + M %*% N)
  return(growth_rates)
}

# wrapper to compute ODE derivatives
Dynamics <- function(time, state, params) {
  dNdt <- state * GetGrowthRates(N = state, pars = params)
  return(list(dNdt))
}

# integrates the dynamics using deSolve
IntegrateDynamics <- function(inistate, pars, endtime, timestep, fn){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars))  
  return(timeseries)
}

# plots time series
PlotSeries <- function(series, title = "HOI Dynamics") {
  meltseries <- melt(series, id.vars = "time")
  pl <- ggplot() +
    geom_line(data = meltseries, aes(x = time, y = value, color = variable),
              size = 1.5, alpha = 0.75) +
    ggtitle(title) +
    xlab("Time") + ylab("Abundance") +
    theme_bw() + theme(legend.position = "none")
  return(pl)
}
