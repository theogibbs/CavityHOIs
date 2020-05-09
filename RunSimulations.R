## Requirements

source("Functions.R")

## Functions

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
  
  ret <- data.frame(SpeciesID = 1:S, Abundances = end_state, Equilibrium = equilibrium, Uninvadeable = uninvadeable)
  print("GetAbd ran! Noice.")
  return(ret)
}

# creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverParams <- function(iterated_params, settings) {
  
  print(paste("Total Rows to Iterate Over:", nrow(iterated_params)))
  
  iterated_params$CommunityID <- 1:nrow(iterated_params)
  pars_list <- alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetAbds, settings)
  
  return(ret_abds)
}

## Parameter Choices

# some scratch space and code to generate the desired combination of parameters
input_S <- 15
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0#seq(0.5, 1, length.out = 2)
input_rho_A <- 0
input_mu_B <- -2
input_sigma_B <- seq(0.5, 0.75, length.out = 2)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

# choosing some of the basic assumptions of the dynamics
settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

# choosing the number of replicates for each parameter combination
num_replicates <- 10
iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))

# setting the number of parallel runs and choosing the current run
arr_length <- 5 # has to agree with the job.slurm file array IDs
cur_ind <- commandArgs(trailingOnly = TRUE)
cur_params <- iterated_params %>%
  mutate(Index = rep(1:arr_length, length.out = nrow(iterated_params))) %>%
  filter(Index == cur_ind)

# running the simulations
system.time(out_abds  <- IterateOverParams(iterated_params = cur_params, settings = settings))

# writing out the data
filename <- "sim505_OnlyPairwise300"
cur_file <- paste0("simdata/", filename, "_", toString(cur_ind), ".csv")
write.csv(out_abds, file = cur_file, row.names = FALSE)

