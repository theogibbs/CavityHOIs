## Packages

library(deSolve)
library(MASS)
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
<<<<<<< HEAD
  scaled_mu_B <- mu_B / S^2
  scaled_sigma_B <- sigma_B / S
=======
  scaled_mu_B <- mu_B / (S^2 - 1)
  scaled_sigma_B <- sigma_B / sqrt(S^2 - 1)
>>>>>>> 7166e0b55f25ddf24a08e541623a207156030fb4
  
  B <- array(rnorm(S*S*S, mean = scaled_mu_B, sd = scaled_sigma_B), c(S, S, S))
  B <- NoSelfHOIs(B)
  return(B)
}

# creates a list of parameters for integrating the dynamics
BuildPars <- function(S = 50, mu_r = 1, sigma_r = 0,mu_d = 1, sigma_d = 0,
                      mu_A = -2, sigma_A = 0.5, rho_A = 0,
                      mu_B = -2, sigma_B = 0.5, rho_B = 0) {
  
  r <- rnorm(S, mean = mu_r, sd = sigma_r)
  d <- rnorm(S, mean = mu_d, sigma_d)
  A <- BuildA(S, mu_A, sigma_A, rho_A)
  B <- BuildB(mu_B, sigma_B, rho_B, mu_A, A)
  
  pars <- list(S = S, r = r, d = d, A = A, B = B)
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
    scale_y_log10() +
    theme_bw() + theme(legend.position = "none")
  return(pl)
}
