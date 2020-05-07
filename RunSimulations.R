## Requirements

source("Functions.R")

##

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
  return(ret)
}

# creates pars and runs the dynamics over the input parameters for a given number of replicates
IterateOverParams <- function(input_params, settings, num_replicates) {
  
  iterated_params <- bind_rows(replicate(num_replicates, input_params, simplify = FALSE))
  iterated_params$CommunityID <- 1:nrow(iterated_params)
  pars_list <- alply(.data = iterated_params, .margins = 1, .fun = BuildPars)
  ret_abds <- ldply(.data = pars_list, .fun = GetAbds, settings)
  
  return(ret_abds)
}



input_S <- 50
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -2
input_sigma_A <- seq(0.5, 1, length.out = 2)
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0#seq(0.5, 1, length.out = 10)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.001, endtime = 1e7, inimin = 0, inimax = 1)

system.time(out_abds <- IterateOverParams(input_params = input_params, settings = settings, num_replicates = 50))


GetStatistics <- function(abds) {
  
  error_bars <- abds %>% group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, CommunityID) %>%
    summarise(CommPhi = unique(sum(Abundances > 0) / S), CommSecAbd = mean(Abundances^2)) %>%
    group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB) %>%
    summarise(Phi = mean(CommPhi), SecAbd = mean(CommSecAbd), ErrorPhi = sd(CommPhi), ErrorSec = sd(CommSecAbd))
  
  abd_stats <- abds %>%
    group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB) %>%
    summarise(Mu = unique(MuA + MuB), Sigma = unique(sqrt(SigmaA^2 + SigmaB^2)),
              MeanAbd = mean(Abundances), FourthAbd = mean(Abundances^4))
  
  ret_stats <- merge(abd_stats, error_bars)
  
  ret_stats <- ret_stats %>%
    mutate(InteractionType = ifelse(SigmaA != 0, ifelse(SigmaB == 0,
                                                        "Pairwise Interactions",
                                                        "Mixed Interactions"),
                                    "Higher Order Interactions")) %>%
    mutate(PlotType = factor(InteractionType,
                             levels = c("Pairwise Interactions", "Mixed Interactions", "Higher Order Interactions")))
  
  
  return(ret_stats)
}

PlotCoexistence <- function(abd_stats) {
  plCoexist <- ggplot(abd_stats, aes(x = Sigma, y = Phi, color = as.factor(Mu))) +
    geom_errorbar(aes(ymin = Phi - 2 * ErrorPhi, ymax = Phi + 2 * ErrorPhi), width = 0) +
    geom_point(size = 2) + theme_bw() + facet_wrap(~ PlotType)
  
  return(plCoexist)
}

PlotCoexistence(GetStatistics(out_abds))
# hard code in a loop over the number of replicates



parallel <- FALSE

if(parallel) {
  arr_length <- 5
  filename <- "sim505_OnlyPairwise300"
  arr_ind <- 1:arr_length
  
  for(cur_ind in arr_ind) {
    #cur_input_params <- input_params %>% filter(ArrInd == cur_ind)
    out_abds <- IterateOverParams(cur_input_params)
    cur_file <- paste0("simdata/", filename, "_", toString(cur_ind), ".csv")
    print(cur_file)
    #write.csv(out_abds, file = cur_file, row.names = FALSE)
  }
  
  
}

#print(system.time(out_data  <- iterate_over_params(input_params, num_replicates = 10)))
#print(head(out_data))

#filename <- "sim505_OnlyPairwise300.csv"
#filename <- paste0("simdata/", filename)
#write.csv(out_data, file = filename, row.names = FALSE)



