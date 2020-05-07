## Requirements

source("Functions.R")

##

# runs the dynamics, gets the endpoints and checks if the system is at equilibrium and is uninvadeable

###### HMMM you could check those things later

GetAbds <- function(pars, inimin = 0, inimax = 1, endtime = 1e7,
                    abd_cutoff = 1e-14, gr_cutoff = 0.01) {
  S <- pars$S
  inistate <- runif(S, min = inimin, max = inimax)
  out <- IntegrateDynamics(inistate, pars, endtime, endtime, Dynamics)
  
  end_state <- as.numeric(out[(nrow(out)),-1])
  survivors <- which(end_state > abd_cutoff)
  invaders <- setdiff(1:S, survivors)
  end_state[invaders] <- 0
  
  growth_rates <- GetGrowthRates(end_state, pars)
  survivor_grs <- growth_rates[survivors]
  equilibrium <- as.logical(prod(survivor_grs <= gr_cutoff))
  
  invader_grs <- growth_rates[invaders]
  pos_invaders <- which(invader_grs > 0)
  uninvadeable <- is_empty(pos_invaders)
  
  ret <- data.frame(SpeciesID = 1:S, Abundances = end_state, Equilibrium = equilibrium, Uninvadeable = uninvadeable)
  return(ret)
}

IterateOverParams <- function(input_params) {
  
  
  
}

iterate_over_params <- function(input_params, num_replicates) {
  coexistence_stats <- get_coexistence_stats(input_params[1,], num_replicates)
  for(i in 2:nrow(input_params)) {
    print(paste("Parameter Combination", i, "of", nrow(input_params)))
    coexistence_stats <- rbind(coexistence_stats, get_coexistence_stats(input_params[i,], num_replicates))
  }
  return(coexistence_stats)
}


# hard code in a loop over the number of replicates
# maybe also have R read in input information from an outside text file???

input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -2
input_sigma_A <- seq(0.5, 1, length.out = 10)
input_rho <- 0
input_mu_B <- 0
input_sigma_B <- 0#seq(0.5, 1, length.out = 10)

abd_cutoff <- 1e-14
gr_cutoff <- 0.001

endtime <- 1e7
inimin <- 0
inimax <- 1

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d, SigmaD = input_sigma_d,
                         MuA = input_mu_A, SigmaA = input_sigma_A, Rho = input_rho, MuB = input_mu_B, SigmaB = input_sigma_B,
                         AbdCutoff = abd_cutoff, GrCutoff = gr_cutoff, EndTime = endtime, IniMin = inimin, IniMax = inimax)

print(system.time(out_data  <- iterate_over_params(input_params, num_replicates = 10)))
print(head(out_data))

filename <- "sim505_OnlyPairwise300.csv"
filename <- paste0("simdata/", filename)
write.csv(out_data, file = filename, row.names = FALSE)



get_coexistence_stats <- function(input_params, num_replicates) {
  
  equilibria <- c()
  invasibility <- c()
  NumCoexist <- c()
  MeanAbd <- c()
  VarAbd <- c()
  
  rep_input_params <- input_params
  
  for(i in 1:num_replicates) {
    pars <- with(input_params, BuildPars(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, Rho, MuB, SigmaB))
    abds_list <- with(input_params, get_abds(pars, IniMin, IniMax, EndTime, AbdCutoff, GrCutoff))
    cur_abd <- abds_list$abds
    equilibria <- c(equilibria, abds_list$equilibrium)
    invasibility <- c(invasibility, abds_list$uninvadeable)
    
    cur_coexist <- sum(abds_list$abds > 0)
    NumCoexist <- c(NumCoexist, cur_coexist)
    MeanAbd <- c(MeanAbd, sum(cur_abd) / cur_coexist)
    cur_abd[cur_abd == 0] <- NA
    VarAbd <- c(VarAbd, var(cur_abd, na.rm = TRUE))
    if(i < num_replicates) {
      rep_input_params <- rbind(rep_input_params, input_params)
    }
  }
  
  ret <- dplyr::bind_cols(rep_input_params, data.frame(Equilibrium = equilibria, Uninvasible = invasibility,
                                                       NumCoexist = NumCoexist, MeanAbd = MeanAbd, VarAbd = VarAbd))
  
  return(ret)
}
