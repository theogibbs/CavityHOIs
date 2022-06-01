## Requirements

source("Functions.R")

## Parameter Choices

# some scratch space and code to generate the desired combination of parameters
input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- 0.25
input_rho_A <- 0
input_mu_B <- c(-1, -3)
input_sigma_B <- seq(0, 2.5, length.out = 10)
input_rho_B <- 0
input_h <- c(1, 3)

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)


# summary of interaction parameters
print(dim(input_params))
print(head(input_params, 10))

# choosing some of the basic assumptions of the dynamics
settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

# choosing the number of replicates to run for a given set of interaction statistics
num_replicates <- 100
iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))
iterated_params$ParsID <- 1:nrow(iterated_params)

# choosing the number of trials to run for a fixed set of interaction parameters but different initial conditions
num_trials <- 1
iterated_params <- bind_rows(replicate(num_trials, iterated_params, simplify = FALSE))
iterated_params$CommunityID <- 1:nrow(iterated_params)


filename <- "sim1207_Saturating30"
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

