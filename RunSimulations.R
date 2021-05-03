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
  
  pars_list <- alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetAbds, settings)
  
  return(ret_abds)
}

## Parameter Choices

# some scratch space and code to generate the desired combination of parameters
input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0#seq(0.5, 1, length.out = 2)
input_rho_A <- 0
input_mu_B <- -2
input_sigma_B <- seq(0.1, 1, length.out = 3)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

input_mu_A <- input_mu_B
input_sigma_A <- input_sigma_B
input_mu_B <- 0
input_sigma_B <- 0

new_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                       SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                       MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

input_params <- rbind(input_params, new_params)

new_params$MuA <- input_mu_A / 2
new_params$MuB <- new_params$MuA
new_params$SigmaA <- new_params$SigmaA / sqrt(2)
new_params$SigmaB <- new_params$SigmaA

input_params <- rbind(input_params, new_params)

print(dim(input_params))
print(head(input_params, 10))

# choosing some of the basic assumptions of the dynamics
settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

# choosing the number of replicates to run for a given set of interaction statistics
num_replicates <- 5
iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choosing the number of trials to run for a gixed set of interaction parameters but different initial conditions
num_trials <- 10
iterated_params <- bind_rows(replicate(num_trials, iterated_params, simplify = FALSE))
iterated_params$CommunityID <- 1:nrow(iterated_params)


filename <- "sim520_MultAtt"
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

