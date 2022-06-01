## Requirements

source("Functions.R")

## Parameter Choices

# some scratch space and code to generate the desired combination of parameters
input_S <- 50
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -2
input_sigma_A <- 0.5
input_rho_A <- 1
input_mu_B <- -4
input_sigma_B <- seq(0.01, 3, length.out = 5)
input_rho_B <- 0
input_h <- 1

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

if(FALSE) {
input_params$ParsID <- 1
pars <- BuildPars(input_params = input_params)
GetGrowthRates(rep(1, times = input_S), pars)
settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.0001, endtime = 1e7, inimin = 0, inimax = 1)

abds_1 <- GetAbds(pars, settings)

settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e6, inimin = 0, inimax = 1)

abds_2 <- GetAbds(pars, settings)

fp_abds <- abds_2$Abundances

out <- with(settings, IntegrateDynamics(abds_2$Initial, pars, endtime, 100, Dynamics))
PlotSeries(out)
fp_sim <- as.numeric(out[nrow(out),(2:ncol(out))]) 

ini_guess <- runif(input_S)
fixed_point <- nleqslv(abds$Abundances, GetGrowthRates, method = "Broyden", jac = NULL, pars)
fixed_point
fval <- sum(fixed_point$fvec)
term_code <- fixed_point$termcd
fp_soln <- fixed_point$x

fp_sim - fp_abds
}
# summary of interaction parameters
print(dim(input_params))
print(head(input_params, 10))

# choosing some of the basic assumptions of the dynamics
settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

# choosing the number of replicates to run for a given set of interaction statistics
num_replicates <- 10
iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choosing the number of trials to run for a fixed set of interaction parameters but different initial conditions
num_trials <- 10
iterated_params <- bind_rows(replicate(num_trials, iterated_params, simplify = FALSE))
iterated_params$CommunityID <- 1:nrow(iterated_params)


filename <- "sim_SatMultAtt"
print(paste0("We are running the simulation: ", filename))

parallel <- FALSE
if(parallel) {
  
  # setting the number of parallel runs and choosing the current run
  arr_length <- 8 # has to agree with the job.slurm file array IDs
  cur_ind <- commandArgs(trailingOnly = TRUE)
  cur_params <- iterated_params %>%
    mutate(Index = rep(1:arr_length, length.out = nrow(iterated_params))) %>%
    filter(Index == cur_ind)
  
  # running the simulations
  print(system.time(out_abds  <- IterateOverParams(iterated_params = cur_params, settings = settings)))
  
  # writing out the data
  cur_file <- paste0("simdata/", filename, "_", toString(cur_ind), ".csv")
  write.csv(out_abds, file = cur_file, row.names = FALSE)
  
} else {
  
  # running the simulations
  print(system.time(out_abds  <- IterateOverParams(iterated_params = iterated_params, settings = settings)))
  
  # writing out the data
  cur_file <- paste0("simdata/", filename, ".csv")
  write.csv(out_abds, file = cur_file, row.names = FALSE)
  
}

