## Dependencies

library(gridExtra)
library(nleqslv)
library(VGAM)
library(viridis)
library(truncnorm)

source("Functions.R")

## Functions
CavitySoln <- function(x, cur_data) {
  
  phi <- x[1]
  avgN <- x[2]
  secN <- x[3]
  v <- x[4]
  
  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * phi^2 * avgN^2))
  var_norm <- with(cur_data,
                   v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + 2 * phi^2 * avgN * secN * SigmaA^2 * RhoB
                          + (SigmaB^2 + RhoB^2 * SigmaA^2) * phi^2 * secN^2))
  
  
  if(var_norm < 0.001) var_norm <- 0.001
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

LabelAbds <- function(abds) {
  
  ret_abds <- abds %>%
    mutate(Interaction= ifelse(SigmaA != 0, ifelse(SigmaB == 0, "Pairwise", "Mixed"), "Higher Order")) %>%
    mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Mixed", "Higher Order"))) %>%
    mutate(Mu = MuA + MuB, Sigma = signif(sqrt(SigmaA^2 + SigmaB^2 + RhoB^2 * SigmaA^2 + RhoB * SigmaA^2)))
  
  return(ret_abds)
}

GetStatistics <- function(abds) {
  ret_stats <- abds %>% group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, CommunityID) %>%
    summarise(CommPhi = unique(sum(Abundances > 0) / length(Abundances)), CommMean = mean(Abundances[Abundances > 0]),
              CommSecAbd = mean(Abundances[Abundances > 0]^2), CommVar = CommSecAbd - CommMean^2,
              CommFourthAbd = mean(Abundances[Abundances > 0]^4), Interaction = unique(Interaction)) %>%
    group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB) %>%
    summarise(Phi = mean(CommPhi), MeanAbd = mean(CommMean), SecAbd = mean(CommSecAbd), VarAbd = mean(CommVar),
              FourthAbd = mean(CommFourthAbd), ErrorPhi = sd(CommPhi), ErrorMean = sd(CommMean), ErrorSec = sd(CommSecAbd),
              ErrorVar = sd(CommVar), ErrorFourth = sd(CommFourthAbd), Mu = unique(MuA + MuB),
              Sigma = unique(signif(sqrt(SigmaA^2 + SigmaB^2))),
              Interaction = unique(Interaction))
  return(ret_stats)
}

PlotCoexistence <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    dplyr::select(Mu, Sigma, RhoA, Phi, MeanAbd, VarAbd, Interaction) %>%
    melt(id.vars = c("Mu", "Sigma", "Interaction", "RhoA"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                           rep("SAD Mean", times = nrow(pred_stats)),
                           rep("SAD Variance", times = nrow(pred_stats)))
  #melt_stats$RhoA <- factor(melt_stats$RhoA, levels = c(-0.5, 0, 0.5),
  #                          ordered = TRUE, labels=c("Negatively Correlated", "Uncorrelated", "Positively Correlated"))
  melt_stats$RhoA <- as.factor(melt_stats$RhoA)
  plCoexist <- ggplot(melt_stats, aes(x = Sigma, y = value, color = RhoA, shape = RhoA)) +
    geom_errorbar(aes(ymin = value - 2 * Error, ymax = value + 2 * Error), width = 0, size = 1) +
    geom_point(size = 4) + theme_bw() + theme(axis.title.y = element_blank(),
                                              legend.position = "top",
                                              text = element_text(size=20),
                                              strip.text.x = element_text(size = 15),
                                              strip.text.y = element_text(size = 15),
                                              legend.text=element_text(size=15)) +
    geom_line(aes(x = Sigma, y = Prediction, color = RhoA), size = 2, alpha = 0.75) +
    facet_wrap(variable ~ Interaction, scales = "free") +
    labs(x = expression("Interaction Heterogeneity"~(sigma)), color = expression("Correlation"~(rho[A])), shape = expression("Correlation"~(rho[A])))
  
  return(plCoexist)
}


PlotCoexistence <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    dplyr::select(Mu, Sigma, RhoA, Phi, MeanAbd, VarAbd, Interaction) %>%
    melt(id.vars = c("Mu", "Sigma", "Interaction", "RhoA"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                           rep("SAD Mean", times = nrow(pred_stats)),
                           rep("SAD Variance", times = nrow(pred_stats)))
  melt_stats$RhoA <- factor(melt_stats$RhoA, levels = c(-0.5, 0, 0.5),
                            ordered = TRUE, labels=c("Negatively Correlated", "Uncorrelated", "Positively Correlated"))
  melt_stats$Mu <- as.factor(melt_stats$Mu)
  melt_stats <- melt_stats %>% filter(variable == "Coexisting Fraction", RhoA == "Uncorrelated")
  
  plCoexist <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
    geom_errorbar(aes(ymin = value - 2 * Error, ymax = value + 2 * Error), width = 0, size = 2) +
    geom_point(size = 7) + theme_bw() + theme(legend.position = "right",
                                              text = element_text(size=30),
                                              strip.text.x = element_text(size = 25),
                                              strip.text.y = element_text(size = 25),
                                              legend.text=element_text(size = 25)) +
    geom_line(aes(x = Sigma, y = Prediction, color = Mu), size = 3, alpha = 0.6) +
    facet_wrap(~ Interaction) + labs(x = expression("Interaction Heterogeneity"),
                                     y = "Coexisting Fraction",
                                     color = "Interaction\nStrength",
                                     shape = "Interaction\nStrength")
  
  return(plCoexist)
}




ParentNorm <- function(x, pred_mean, pred_sec) {
  mu <- x[1]
  sigma <- x[2]
  pred_var <- pred_sec - pred_mean^2
  
  mu_eq <- pred_mean - etruncnorm(a = 0, b = Inf, mean = mu, sd = sigma)
  var_eq <- pred_var - vtruncnorm(a = 0, b = Inf, mean = mu, sd = sigma)
  
  return(c(mu_eq, var_eq))
}

PlotHist <- function(proc_abds, hist_sigmas) {
  plot_abds <- proc_abds %>%
    filter(Sigma == hist_sigmas[1] | Sigma == hist_sigmas[2]) %>%
    filter(Abundances > 0)
  
  pred_data <- plot_abds %>%
    dplyr::select(-one_of("Abundances", "SpeciesID", "CommunityID", "Index")) %>%
    unique()
  
  pred_abds <- GetPredictions(pred_data)
  pred_abds$AdjPredMean <- 0
  pred_abds$PredSD <- 0

  for(i in 1:nrow(pred_abds)) {
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    new_preds <- nleqslv(c(pred_mean, sqrt(pred_sec - pred_mean^2)), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$PredSD[i] <- new_preds[2]
  }
  
  plot_abds <- merge(plot_abds, pred_abds)
  plot_abds <- plot_abds %>%
    mutate(PredDensity = dtruncnorm(Abundances, a = 0, b = Inf, mean = AdjPredMean, sd = PredSD))
  
  plHist <- ggplot(plot_abds, aes(x = Abundances, y = ..density..)) +
    geom_histogram(fill = "white", color = "black", binwidth = 0.02) +
    #ggtitle("Species Abundance Distributions") +
    facet_grid(Sigma ~ Interaction, labeller = label_bquote(rows = sigma == .(Sigma))) +
    theme_bw() + theme(strip.text.x = element_text(size = 25),
                       strip.text.y = element_text(size = 25),
                       text = element_text(size=30)) + ylab("Density") +
    geom_line(aes(x = Abundances, y = PredDensity), color = "blue", size = 1)
  return(plHist)
}

PlotHist <- function(proc_abds, hist_sigmas) {
  plot_abds <- proc_abds %>%
    filter(Sigma == hist_sigmas[1] | Sigma == hist_sigmas[2]) %>%
    filter(Abundances > 0)
  
  pred_data <- plot_abds %>%
    dplyr::select(-one_of("Abundances", "SpeciesID", "CommunityID", "Index")) %>%
    unique()
  
  pred_abds <- GetPredictions(pred_data)
  pred_abds$AdjPredMean <- 0
  pred_abds$PredSD <- 0
  
  for(i in 1:nrow(pred_abds)) {
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    new_preds <- nleqslv(c(pred_mean, sqrt(pred_sec - pred_mean^2)), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$PredSD[i] <- new_preds[2]
  }
  
  plot_abds <- merge(plot_abds, pred_abds)
  plot_abds <- plot_abds %>%
    mutate(PredDensity = dtruncnorm(Abundances, a = 0, b = Inf, mean = AdjPredMean, sd = PredSD))
  
  plHist <- ggplot(plot_abds, aes(x = Abundances, y = ..density..)) +
    geom_histogram(fill = "white", color = "black", binwidth = 0.02) +
    #ggtitle("Species Abundance Distributions") +
    facet_grid(Mu ~ Sigma, labeller = label_bquote(rows = mu[B] == .(Mu), cols = sigma[B] == .(Sigma))) +
    theme_bw() + theme(strip.text.x = element_text(size = 25),
                       strip.text.y = element_text(size = 25),
                       text = element_text(size=30)) + ylab("Density") +
    geom_line(aes(x = Abundances, y = PredDensity), color = "blue", size = 1)
  return(plHist)
}

fileregex <- "sim511_CorrHOIs."
corr_abds <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))


fileregex <- "sim508_OnlyPairwise300."
out_abds_pw <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

fileregex <- "sim510_OnlyHOIs300."
out_abds_hois <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

fileregex <- "sim511_Mixed300."
out_abds_mx <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

out_abds_og <- rbind(out_abds_pw, out_abds_hois, out_abds_mx)


fileregex <- "sim51._AllInts300."
out_abds_all <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))


fileregex <- "sim522_PairCorr."
out_abds_pair <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

out_abds_og <- rbind(out_abds_og, out_abds_all)
out_abds_og$ParsID <- 0
out_abds <- rbind(out_abds_og, out_abds_pair)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
#proc_abds <- LabelAbds(out_abds) %>% filter(MuB == 0 & SigmaB == 0 | MuA == 0 & SigmaA == 0)
proc_abds <- LabelAbds(out_abds) # %>% filter(Mu == -4)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)
#for(i in 1:nrow(pred_stats)) {
#  cur_mean <- pred_stats$PredMean[i]
#  cur_sec <- pred_stats$PredSec[i]
#  new_preds <- nleqslv(c(cur_mean, sqrt(cur_sec - cur_mean^2)), fn = ParentNorm, jac = NULL, cur_mean, cur_sec)$x
#  pred_stats$AdjPredMean[i] <- new_preds[1]
#  pred_stats$PredSD[i] <- new_preds[2]
#}


plCoexist <- PlotCoexistence(pred_stats)
show(plCoexist)


#png("../CavityHOIs-Notes/NegcorrCoexistence.png", width = 3500, height = 1500, res = 300)
#plCoexist
#dev.off()

png("../CavityHOIs-Notes/Coexistence.png", width = 3500, height = 3000, res = 300)
plCoexist
dev.off()

plHist <- PlotHist(LabelAbds(out_abds_og) %>% filter(Interaction %in% c("Pairwise", "Higher Order")), c(0.5, 1))
show(plHist)
png("../CavityHOIs-Notes/Histogram.png", width = 2000, height = 1500, res = 300)
plHist
dev.off()

## Heatmap plotting

GetPredictions <- function(out_data, max_trials) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredSec = c(), PredV = c(), FuncVal = c(), TermCode = c())
  
  for(row in 1:nrow(out_data)) {
    
    if(row %% 10 == 0) {
      print("Row number:")
      print(row)
    }
    cur_data <- out_data[row,]
    ini_guess <- runif(4, 0, 1)
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


# some scratch space and code to generate the desired combination of parameters

input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0#seq(-2, -4, length.out = 50)
input_sigma_A <- 0#seq(0.5, 1, length.out = 2)
input_rho_A <- c(-0.5, 0, 0.5)
input_mu_B <- seq(0, -5, length.out = 20)
input_sigma_B <- seq(0.01, 1.5, length.out = 20)
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

new_params$MuA <- new_params$MuA / 2
new_params$MuB <- new_params$MuA
new_params$SigmaA <- new_params$SigmaA / sqrt(2)
new_params$SigmaB <- new_params$SigmaA

input_params <- rbind(input_params, new_params)

plot_pred <- GetPredictions(LabelAbds(input_params), 10) %>% mutate(PredFraction = ifelse(TermCode == 1, PredFraction, NA))

plHeatMap <- ggplot(plot_pred, aes(x = Mu, y = Sigma, fill = PredFraction)) +
  geom_tile() + scale_fill_viridis(end = 0.9, na.value="darkred") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme(strip.background = element_rect(fill="lightgray", color = "black"),
                                               aspect.ratio = 1,
                                               axis.text = element_text( size = 10 ),
                                               panel.background = element_rect(fill = NA),
                                               text = element_text(size=20),
                                               legend.key.size = unit(3, "cm"),
                                               legend.key.width = unit(1,"cm")) +
  facet_grid(RhoA~as.factor(Interaction), labeller = label_bquote(rows = rho[A] == .(RhoA))) +
  coord_fixed(ratio = (diff(range(plot_pred$Mu)) / diff(range(plot_pred$Sigma)))) +
  labs(x = expression("Mean Interaction Strength"~(mu)), y = expression("Interaction Heterogeneity"~(sigma)),
       fill = expression("Coexisting\nFraction"*(phi)))
show(plHeatMap)
png("../CavityHOIs-Notes/HeatMap.png", width = 3000, height = 3000, res = 300)
plHeatMap
dev.off()

ratio_preds <- plot_pred %>% select(Mu, Sigma, RhoA, PredFraction, Interaction)

pw_preds <- ratio_preds %>% filter(Interaction == "Pairwise")
hoi_preds <- ratio_preds %>% filter(Interaction == "Higher Order")

ratio_preds <- merge(pw_preds, hoi_preds, by = c("Mu", "Sigma", "RhoA")) %>%
  mutate(Ratio = PredFraction.x/PredFraction.y)
 # mutate(Ratio = ifelse((!is.na(PredFraction.x) & is.na(PredFraction.y)), max(PredFraction.x/PredFraction.y, na.rm = TRUE), PredFraction.x/PredFraction.y))

plRatio <- ggplot(ratio_preds, aes(x = Mu, y = Sigma, fill = Ratio)) +
  geom_tile() + scale_fill_gradient2(midpoint = 1, na.value = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme(strip.background = element_rect(fill="lightgray", color = "black"),
                                               aspect.ratio = 1,
                                               axis.text = element_text( size = 10 ),
                                               panel.background = element_rect(fill = NA),
                                               text = element_text(size=20),
                                               legend.key.size = unit(3, "cm"),
                                               legend.key.width = unit(1,"cm")) +
  facet_grid(~RhoA, labeller = label_bquote(rows = rho[yA] == .(RhoA))) +
  coord_fixed(ratio = (diff(range(plot_pred$Mu)) / diff(range(plot_pred$Sigma)))) +
  labs(x = expression("Mean Interaction Strength"~(mu)), y = expression("Interaction Heterogeneity"~(sigma)),
       fill = expression("Ratio"))
plRatio

## correlated communities

PlotCorr <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    select(Mu, Sigma, Phi, MeanAbd, VarAbd, RhoB, Interaction) %>%
    melt(id.vars = c("Mu", "Sigma", "Interaction", "RhoB"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)), rep("SAD Mean", times = nrow(pred_stats)),
                           rep("SAD Variance", times = nrow(pred_stats)))
  
  plCoexist <- ggplot(melt_stats, aes(x = Sigma, y = value, color = as.factor(RhoB), shape = as.factor(RhoB))) +
    geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 1) +
    geom_point(size = 4) + theme_bw() + theme(axis.title.y=element_blank(),
                                              text = element_text(size=15),
                                              strip.text.x = element_text(size = 15),
                                              strip.text.y = element_text(size = 15),
                                              legend.text=element_text(size=10)) + labs(color = "RhoB", shape = "RhoB") + 
    geom_line(aes(x = Sigma, y = Prediction, color = as.factor(RhoB)), size = 2, alpha = 0.75) +
    facet_wrap(~ variable, scales = "free")
  
  return(plCoexist)
}

CavitySoln <- function(x, cur_data) {
  
  phi <- x[1]
  avgN <- x[2]
  secN <- x[3]
  v <- x[4]
  
  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * phi^2 * avgN^2))
  var_norm <- with(cur_data,
                   v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + 2 * phi^2 * avgN * secN * SigmaA^2 * RhoB
                          + (SigmaB^2 + RhoB^2 * SigmaA^2) * phi^2 * secN^2 + 1 *phi^3 * SigmaA^2 * RhoB^2 * secN * avgN^2))
  
  
  if(var_norm < 0.001) var_norm <- 0.001
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


fileregex <- "sim51._CorrHOIs."
out_abds_corr <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

fileregex <- "sim511_Mixed300."
out_abds_mx <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

corr_abds <- rbind(out_abds_corr, out_abds_mx)

CheckAbds(corr_abds)
corr_abds <- corr_abds %>% filter(Equilibrium == TRUE)

proc_abds_corr <- LabelAbds(corr_abds)

abd_stats <- GetStatistics(proc_abds_corr)
pred_stats <- GetPredictions(abd_stats, 100)

pred_stats <- pred_stats %>% mutate(Phi = ifelse(RhoB == 0, Phi, Phi / 2)) %>% mutate(Phi = ifelse(Phi < 0.5, 2 * Phi, Phi))


plCorr <- PlotCorr(pred_stats)
show(plCorr)
png("../CavityHOIs-Notes/Correlation.png", width = 3000, height = 1000, res = 300)
plCorr
dev.off()


## Multiple Attractor Analysis

fileregex <- "sim520_MultAtt_."
out_abds_ma <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

CheckAbds(out_abds_ma)
proc_abds_ma <- LabelAbds(out_abds_ma)

ma_data <- proc_abds_ma %>% group_by(SpeciesID, S, Mu, Sigma, Interaction, ParsID) %>%
  summarise(UFP = isTRUE(all.equal(min(Abundances), max(Abundances), tolerance = 1e-5))) %>%
  group_by(S, Mu, Sigma, Interaction) %>% summarise(UFP = prod(UFP))
print(ma_data)

## getting the non-phi adjusted mean and variance

GetAdjPreds <- function(pred_data) {
  pred_abds <- pred_data
  pred_abds$AdjPredMean <- 0
  pred_abds$AdjPredSD <- 0
  pred_abds$PredSD <- sqrt(pred_abds$PredSec - pred_abds$PredMean^2)
  
  for(i in 1:nrow(pred_abds)) {
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    new_preds <- nleqslv(c(pred_mean, sqrt(pred_sec - pred_mean^2)), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$AdjPredSD[i] <- new_preds[2]
  }
  
  pred_abds$ExpMin <- pred_abds$PredMean - pred_abds$PredSD * sqrt(2 * log(pred_abds$S))
  pred_abds$AdjExpMin <- pred_abds$AdjPredMean - pred_abds$AdjPredSD * sqrt(2 * log(pred_abds$S))
  
  return(pred_abds)
}

adj_pred_stats <- GetAdjPreds(pred_stats)

melt_adj <- adj_pred_stats %>%
  ungroup() %>%
  select(Sigma, RhoA, PredMean, PredSD, AdjPredSD, AdjPredMean, ExpMin, AdjExpMin, Interaction) %>%
  melt(id.vars = c("Sigma", "Interaction", "RhoA")) %>%
  mutate(type = ifelse(variable %in% c("PredMean", "AdjPredMean"), "Mean", ifelse(variable %in% c("PredSD", "AdjPredSD"), "Standard Deviation", "Expected Minimum"))) %>%
  mutate(var = ifelse(variable %in% c("PredMean", "PredSD", "ExpMin"), "Predicted", "Adjusted"))

plAdj <- ggplot(melt_adj, aes(x = Sigma, y = value, color = Interaction, linetype = var)) +
  geom_line(size = 2) + theme_bw() + theme(axis.title.y = element_blank(),
                                           legend.position = "top",
                                           text = element_text(size=20),
                                           strip.text.x = element_text(size = 15),
                                           strip.text.y = element_text(size = 15),
                                           legend.text=element_text(size=15)) +
  facet_wrap(type ~ RhoA, scales = "free")
plAdj

png("../CavityHOIs-Notes/Adjusted.png", width = 3000, height = 3000, res = 300)
plAdj
dev.off()

## rescaling

# some scratch space and code to generate the desired combination of parameters

input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- seq(-2, -4, length.out = 15)
input_sigma_A <- seq(0, 1, length.out = 50)
input_rho_A <- 0
input_mu_B <- 0#seq(-2, -5, length.out = 50)
input_sigma_B <- 0#seq(0.01, 1, length.out = 50)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)
input_params$InputMu <- input_params$MuA
input_params$InputSigma <- input_params$SigmaA


plot_pred <- GetPredictions(LabelAbds(input_params), 100)

pw_pred <- plot_pred$PredFraction
pw_mean <- plot_pred$PredMean
pw_sec <- plot_pred$PredSec

input_mu <- input_params$InputMu
input_sigma <- input_params$InputSigma
input_mu_B <- input_mu * (input_mu_d - input_mu) / input_mu_r
input_mu_B <- input_mu / pw_mean / pw_pred
input_sigma_B <- sqrt(input_sigma^2 * (1 - input_sigma^2) * (input_mu_d - input_mu)^2 / input_mu_r^2)
input_sigma_B <- sqrt(input_sigma^2 / pw_sec / pw_pred)
input_mu_A <- 0
input_sigma_A <- 0

new_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                       SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                       MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

new_params <- as.data.frame(cbind(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                       SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                       MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B))
new_params$InputMu <- input_mu
new_params$InputSigma <- input_sigma

input_params <- rbind(input_params, new_params)

plot_pred <- GetPredictions(LabelAbds(input_params), 100) %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
 # mutate(Mu = ifelse(Interaction == "Pairwise", PredFraction * Mu, Mu * (-1 + 2 * PredFraction))) %>%
  mutate(Mu = as.factor(Mu))

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(InputMu, InputSigma, RhoA, Interaction, PredFraction, PredMean, PredVar) %>%
  melt(id.vars = c("InputMu", "InputSigma", "Interaction", "RhoA"))

plPred <- ggplot(melt_pred, aes(x = InputSigma, y = value, color = as.factor(InputMu), shape = Interaction)) +
  geom_point(size = 4, alpha = 0.5) + theme_bw() + theme(axis.title.y = element_blank(),
                                           legend.position = "top",
                                           text = element_text(size=20),
                                           strip.text.x = element_text(size = 15),
                                           strip.text.y = element_text(size = 15),
                                           legend.text=element_text(size=15)) +
  facet_wrap(~ variable, scales = "free")
plPred

plMu <-  ggplot(input_params %>% filter(SigmaB != 0), aes(x = InputMu, y = MuB, color = InputMu)) +
  geom_point(size = 4) + theme_bw() + geom_abline(slope = 1, intercept = 0)
plMu
plSigma <- ggplot(input_params %>% filter(SigmaB != 0), aes(x = InputSigma, y = SigmaB, color = InputMu)) +
  geom_point(size = 4) + theme_bw() + geom_abline(slope = 1, intercept = 0)
plSigma



########## DRAFTING FIGURE 2


GetAdjPreds <- function(pred_data) {
  pred_abds <- pred_data
  pred_abds$AdjPredMean <- 0
  pred_abds$AdjPredSD <- 0
  pred_abds$PredSD <- sqrt(pred_abds$PredSec - pred_abds$PredMean^2)
  
  for(i in 1:nrow(pred_abds)) {
    pred_mean <- pred_abds$PredMean[i]
    pred_sec <- pred_abds$PredSec[i]
    new_preds <- nleqslv(c(pred_mean, sqrt(pred_sec - pred_mean^2)), fn = ParentNorm, jac = NULL, pred_mean, pred_sec)$x
    pred_abds$AdjPredMean[i] <- new_preds[1]
    pred_abds$AdjPredSD[i] <- new_preds[2]
  }
  
  pred_abds$ExpMin <- pred_abds$PredMean - pred_abds$PredSD * sqrt(2 * log(pred_abds$S))
  pred_abds$AdjExpMin <- pred_abds$AdjPredMean - pred_abds$AdjPredSD * sqrt(2 * log(pred_abds$S))
  
  return(pred_abds)
}


input_S <- 50
input_mu_r <- 10
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0
input_rho_A <- 0
input_mu_B <- c(-0.5, -1, -2) * 10
input_sigma_B <- seq(1e-7, -max(input_mu_B) / sqrt(3) / input_S, length.out = 50)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

plot_pred <- GetPredictions(LabelAbds(input_params), 10) %>% filter(TermCode == 1)
plot_pred <- GetAdjPreds(plot_pred)

plot_pred <- plot_pred %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
  mutate(Mu = as.factor(Mu)) %>%
  mutate(AdjRatio = AdjPredMean / sqrt(2) / AdjPredSD) %>%
  group_by(SigmaR) %>%
  mutate(AdjRatio =  AdjRatio / max(AdjRatio))

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(MuB, SigmaB, RhoA, Interaction, SigmaR, PredFraction, AdjPredMean, AdjPredSD, AdjExpMin, AdjRatio) %>%
  mutate(AdjPredSD = AdjPredSD^2) %>%
  melt(id.vars = c("MuB", "SigmaB", "Interaction", "RhoA", "SigmaR")) %>%
  filter(variable != "AdjExpMin") %>%
  mutate(variable = factor(variable, levels = c("PredFraction", "AdjPredMean", "AdjPredSD", "AdjRatio"), labels = c("Coexisting\nFraction",
                                                                                                        "Adjusted\nMean", "Adjusted\nVariance", "Mean to SD\nRatio")))

plPred <- ggplot(melt_pred, aes(x = SigmaB, y = value, color = as.factor(MuB))) +
  geom_line(size = 2, alpha = 0.75) + theme_bw() + theme(axis.title.y = element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.placement = "outside",
                                                         text = element_text(size=25),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size=25),
                                                         panel.spacing = unit(2, "lines")) +
  facet_grid(variable ~ SigmaR, labeller = label_bquote(cols = sigma[R] == .(SigmaR)),
             scales = "free", switch = "y") +
  labs(x = expression("HOI Heterogeneity"~(sigma[B])), color = expression("HOI\nStrength"~(mu[B]))) +
  ggtitle("Coexistence for Higher Order Interactions")
plPred

jpeg("../CavityHOIs-Paper/Fig3.jpeg", width = 4500, height = 3500, res = 300)
plPred
dev.off()

### purely pairwise analog


input_S <- 300
input_mu_r <- 1
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- c(0, -0.5, -1)
input_sigma_A <- seq(0.01, 0.75, length.out = 50)
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

plot_pred <- GetPredictions(LabelAbds(input_params), 100) %>% filter(TermCode == 1)
plot_pred <- GetAdjPreds(plot_pred)

plot_pred <- plot_pred %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
  mutate(Mu = as.factor(Mu)) %>%
  mutate(AdjRatio = AdjPredMean / sqrt(2) / AdjPredSD) %>%
  group_by(SigmaR) %>%
  mutate(AdjRatio =  AdjRatio / max(AdjRatio))

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(MuA, SigmaA, RhoA, Interaction, SigmaR, PredFraction, AdjPredMean, AdjPredSD, AdjExpMin, AdjRatio) %>%
  mutate(AdjPredSD = AdjPredSD^2) %>%
  melt(id.vars = c("MuA", "SigmaA", "Interaction", "RhoA", "SigmaR")) %>%
  filter(variable != "AdjExpMin") %>%
  mutate(variable = factor(variable, levels = c("PredFraction", "AdjPredMean", "AdjPredSD", "AdjRatio"), labels = c("Coexisting\nFraction",
                                                                                                                    "Adjusted\nMean", "Adjusted\nVariance",
                                                                                                                    "Mean to SD\nRatio")))

plPred <- ggplot(melt_pred, aes(x = SigmaA, y = value, color = as.factor(MuA))) +
  geom_line(size = 2, alpha = 0.75) + theme_bw() + theme(axis.title.y = element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.placement = "outside",
                                                         text = element_text(size=25),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size=25),
                                                         panel.spacing = unit(2, "lines")) +
  facet_grid(variable ~ SigmaR, labeller = label_bquote(cols = sigma[R] == .(SigmaR)),
             scales = "free", switch = "y") +
  labs(x = expression("Interaction Heterogeneity"~(sigma[A])), color = expression("Interaction\nStrength"~(mu[A]))) +
  ggtitle("Coexistence for Pairwise Interactions")
plPred


jpeg("../CavityHOIs-Paper/Fig4.jpeg", width = 4500, height = 3500, res = 300)
plPred
dev.off()

### varying mean interaction strength for pairwise and higher order parameterizations

input_S <- 300
input_mu_r <- 1
input_sigma_r <- c(0, 0.25)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- seq(0, -5, length.out = 50)
input_sigma_A <- 0.5
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)


input_mu_B <- input_mu_A
input_mu_A <- 0
input_sigma_B <- input_sigma_A
input_sigma_A <- 0

input_params <- rbind(input_params, crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B))

plot_pred <- GetPredictions(LabelAbds(input_params), 100) %>% filter(TermCode == 1)
plot_pred <- GetAdjPreds(plot_pred)

plot_pred <- plot_pred %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
  mutate(AdjRatio = AdjPredMean / sqrt(2) / AdjPredSD)

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(Mu, Sigma, Interaction, SigmaR, PredFraction, AdjPredMean, AdjPredSD, AdjExpMin, AdjRatio) %>%
  mutate(AdjPredSD = AdjPredSD^2) %>%
  melt(id.vars = c("Mu", "Sigma", "Interaction", "SigmaR")) %>%
  filter(variable != "AdjExpMin") %>%
  mutate(Interaction = ifelse(Interaction == "Higher Order", "Higher\nOrder", "Pairwise")) %>%
  mutate(variable = factor(variable, levels = c("PredFraction", "AdjPredMean", "AdjPredSD", "AdjRatio"), labels = c("Coexisting\nFraction",
                                                                                                                    "Adjusted\nMean", "Adjusted\nVariance",
                                                                                                                    "Mean to SD\nRatio")))

plPred <- ggplot(melt_pred, aes(x = Mu, y = value, color = Interaction)) +
  geom_line(size = 2, alpha = 0.75) + theme_bw() + theme(axis.title.y = element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.placement = "outside",
                                                         text = element_text(size=25),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size=25),
                                                         panel.spacing = unit(2, "lines")) +
  scale_x_reverse() +
  facet_grid(variable ~ SigmaR, labeller = label_bquote(cols = sigma[R] == .(SigmaR)),
             scales = "free", switch = "y") +
  labs(x = expression("Interaction Strength"~(mu)), color = expression("Interaction\nType")) +
  ggtitle("Coexistence for Varying Interaction Strengths")
plPred


### other stuff

input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- c(-1, -2, -4)
input_sigma_A <- seq(0.5, 0.75, length.out = 5)
input_rho_A <- c(-0.5, 0, 0.5)
input_mu_B <- seq(0, -5, length.out = 20)
input_sigma_B <- c(0.01, 0.5)#seq(0.01, 1, length.out = 25)
input_rho_B <- 0


plot_pred <- GetPredictions(LabelAbds(input_params), 100) %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
  mutate(Mu = as.factor(Mu))

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(MuA, SigmaA, MuB, SigmaB, RhoA, Interaction, PredFraction, PredMean, PredVar) %>%
  melt(id.vars = c("MuA", "SigmaA", "MuB", "SigmaB", "Interaction", "RhoA")) %>%
  filter(variable == "PredFraction")


input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)


ggplot(melt_pred, aes(x = MuB, y = value, shape = as.factor(SigmaA), color = as.factor(MuA))) +
  geom_point(size = 4, alpha = 0.5) + theme_bw() + theme(axis.title.y = element_blank(),
                                                         legend.position = "top",
                                                         text = element_text(size=20),
                                                         strip.text.x = element_text(size = 15),
                                                         strip.text.y = element_text(size = 15),
                                                         legend.text=element_text(size=15)) +
  facet_grid(SigmaB ~ RhoA)


# WORSE ATTEMPTS BELOW I BELIEVE

# some scratch space and code to generate the desired combination of parameters

input_S <- 300
input_mu_r <- seq(4, 6, length.out = 3)
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -4#seq(-2, -4, length.out = 50)
input_sigma_A <- seq(0.1, 0.5, length.out = 20)
input_rho_A <- 0
input_mu_B <- 0#seq(-2, -5, length.out = 50)
input_sigma_B <- 0#seq(0.01, 1, length.out = 50)
input_rho_B <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)


input_mu_B <- input_mu_A
input_sigma_B <- input_sigma_A
input_mu_A <- 0
input_sigma_A <- 0

new_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                       SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                       MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B)

input_params <- rbind(input_params, new_params)

plot_pred <- GetPredictions(LabelAbds(input_params), 100) %>%
  mutate(PredVar = PredSec - PredMean^2) %>%
  # mutate(Mu = ifelse(Interaction == "Pairwise", PredFraction * Mu, Mu * (-1 + 2 * PredFraction))) %>%
  mutate(Mu = as.factor(Mu))

melt_pred <- plot_pred %>%
  ungroup() %>%
  dplyr::select(Mu, MuR, Sigma, RhoA, Interaction, PredFraction, PredMean, PredVar) %>%
  melt(id.vars = c("Mu", "MuR", "Sigma", "Interaction", "RhoA"))

plPred <- ggplot(melt_pred, aes(x = Sigma, y = value, color = as.factor(MuR), shape = Interaction)) +
  geom_point(size = 4, alpha = 0.5) + theme_bw() + theme(axis.title.y = element_blank(),
                                                         legend.position = "top",
                                                         text = element_text(size=20),
                                                         strip.text.x = element_text(size = 15),
                                                         strip.text.y = element_text(size = 15),
                                                         legend.text=element_text(size=15)) +
  facet_wrap(~ variable, scales = "free")
plPred

grid.arrange(plMu, plSigma, nrow = 1)




