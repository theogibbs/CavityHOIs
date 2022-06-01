
source("Functions.R")

## Functions

GetCriticalSigma <- function(pars, settings, sigmas) {
  
  max_trials <- length(sigmas)
  coexist <- TRUE
  cur_trial <- 1
  
  pars$A <- pars$A / ifelse(sd(pars$A) == 0, 1, sd(pars$A))
  pars$B <- pars$B / ifelse(sd(pars$B) == 0, 1, sd(pars$B))
  
  while(coexist && (cur_trial < max_trials)) {
    
    cur_pars <- pars
    sigma <- sigmas[cur_trial]
    cur_pars$A <- sigma * pars$A
    cur_pars$B <- sigma * pars$B

    cur_abds <- GetAbds(cur_pars, settings)
    coexist <- prod(cur_abds$Abundances > 0)
    cur_eq <- unique(cur_abds$Equilibrium)
    cur_uninv <- unique(cur_abds$Uninvadeable)
    cur_trial <- cur_trial + 1
  }
  
  out_sigma <- mean(c(sigmas[cur_trial-2], sigma))
  
  return(data.frame(Sigma = out_sigma, Equilibrium = cur_eq, Uninvadeable = cur_uninv, MaxSteps = (cur_trial == max_trials)))
}

## Parameter Choices

# some scratch space and code to generate the desired combination of parameters
input_S <- seq(5, 40, by = 3)
input_mu_r <- 1.5
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0
input_sigma_B <- 1
input_rho_A <- 0
input_mu_B <- 0
input_rho_B <- 0
input_h <- 0
in_int <- "Higher Order"

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h, Interaction = in_int)

input_sigma_B <- 0
input_sigma_A <- 1
in_int <- "Pairwise"

input_params <- rbind(input_params, crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                                             SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                                             MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h, Interaction = in_int))


# choosing the number of replicates to run for a given set of interaction statistics
num_replicates <- 10
iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choosing the number of trials to run for a fixed set of interaction parameters but different initial conditions
num_trials <- 1
iterated_params <- bind_rows(replicate(num_trials, iterated_params, simplify = FALSE))
iterated_params$CommunityID <- 1:nrow(iterated_params)

# choosing the sensitivity, upper bound and maximum number of trials for the critical sigma search
max_sigma <- 1
num_sigmas <- 2e3
step_size <- max_sigma / num_sigmas
print(step_size)
test_sigmas <- seq(0, max_sigma, length.out = num_sigmas)

settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e6, inimin = 0, inimax = 1)

what_time_is_it <- Sys.time()
out_sigmas <- data.frame()
for(i in 1:nrow(iterated_params)) {
  print(paste("Run", i, "out of", nrow(iterated_params)))
  cur_params <- iterated_params[i,]
  pars <- BuildPars(cur_params)
  cur_sigmas <- GetCriticalSigma(pars, settings, test_sigmas)
  cur_sigmas <- cbind(cur_params, cur_sigmas)
  out_sigmas <- rbind(out_sigmas, cur_sigmas)
}

print(Sys.time() - what_time_is_it)

file_name <- "sim1204_MayBound.csv"

write.csv(out_sigmas, file = paste0("simdata/", file_name), row.names = FALSE)




