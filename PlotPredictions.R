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
  fourthN <- x[4]
  v <- x[5]
  
  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * phi^2 * avgN^2))
  var_norm <- with(cur_data, v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + SigmaB^2* phi^2 * secN^2 ))
  
  
  if(var_norm < 0.001) var_norm <- 0.001
  error_fn <- erf(avg_norm / sqrt(2 * var_norm))
  
  eq1 <- 0.5 * (1 + error_fn)
  eq2 <- avg_norm / 2 * (1 + error_fn) + sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  eq3 <- (avg_norm^2 + var_norm) / 2 * (1 + error_fn)
  eq3 <- eq3 + avg_norm * sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  eq4 <- (3 * var_norm^2 + 6 * avg_norm^2 * var_norm + avg_norm^4) / 2 * (1 + error_fn)
  eq4 <- eq4 + avg_norm * sqrt(var_norm) * (5 * avg_norm * var_norm + avg_norm^3) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  
  eq1 <- phi - eq1
  eq2 <- avgN - (1 / phi) * eq2
  eq3 <- secN - (1 / phi) * eq3
  eq4 <- fourthN - (1 / phi) * eq4
  eq5 <- with(cur_data, 1 - v * (MuD - 0 * phi * MuB * avgN / S + phi * RhoA * SigmaA^2 * v))
  
  return(c(eq1, eq2, eq3, eq4, eq5))
}


GetPredictions <- function(out_data) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredSec = c(), PredFourth = c(), PredV = c())
  
  for(row in 1:nrow(out_data)) {
    
    if(row %% 100 == 0) {
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
    ini_var <- ini_sec - ini_mean^2
    
    ini_fourth <- 3 * ini_var^2 + 6 * ini_mean^2 * ini_var + ini_mean^4
    
    h <- function(ini_v, cur_data) with(cur_data, 1 - ini_v * (MuD - 2 * MuB * ini_mean / (S+1) + RhoA * SigmaA^2 * ini_v))
    for(curx in testxs) hxs <- c(hxs, h(curx, cur_data))
    max_v <- testxs[hxs < 0][1]
    ini_v <- uniroot(h, c(0, max_v), cur_data)$root
    
    ini_guess <- c(1 - 0.25 * with(cur_data, sqrt(SigmaA^2 + SigmaB^2)), ini_mean, ini_sec, ini_fourth, ini_v)
    cav_soln <- nleqslv(ini_guess, CavitySoln, method = "Broyden", jac = NULL, cur_data)$x

    PRED <- rbind(PRED, data.frame(PredFraction = cav_soln[1], PredMean = cav_soln[2],
                                   PredSec = cav_soln[3], PredFourth = cav_soln[4], PredV = cav_soln[5]))
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
    mutate(Mu = MuA + MuB, Sigma = signif(sqrt(SigmaA^2 + SigmaB^2)))
  
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
              Sigma = unique(sqrt(SigmaA^2 + SigmaB^2)), Interaction = unique(Interaction))
  return(ret_stats)
}

PlotCoexistence <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    select(Mu, Sigma, Phi, MeanAbd, VarAbd, Interaction) %>%
    melt(id.vars = c("Mu", "Sigma", "Interaction"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("(a) Coexisting Fraction", times = nrow(pred_stats)), rep("(b) SAD Mean", times = nrow(pred_stats)),
                           rep("(c) SAD Variance", times = nrow(pred_stats)))
  
  plCoexist <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Interaction, shape = Interaction)) +
    geom_errorbar(aes(ymin = value - 2 * Error, ymax = value + 2 * Error), width = 0, size = 1) +
    geom_point(size = 4) + theme_bw() + theme(axis.title.y=element_blank(),
                                              legend.title = element_blank(),
                                              text = element_text(size=15),
                                              strip.text.x = element_text(size = 15),
                                              strip.text.y = element_text(size = 15),
                                              legend.text=element_text(size=10)) +
    geom_line(aes(x = Sigma, y = Prediction, color = Interaction), size = 2, alpha = 0.75) +
    facet_wrap(~ variable, scales = "free") + labs(x = expression("Interaction Heterogeneity"~(sigma)))
    
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
    select(-one_of("Abundances", "SpeciesID", "CommunityID", "Index")) %>%
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
    ggtitle("Species Abundance Distributions") +
    facet_grid(Sigma ~ Interaction, labeller = label_bquote(rows = sigma == .(Sigma))) +
    theme_bw() + theme(strip.text.x = element_text(size = 15),
                       strip.text.y = element_text(size = 15),
                       text = element_text(size=20)) + ylab("Density") +
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

out_abds <- rbind(out_abds_pw, out_abds_hois, out_abds_mx)


fileregex <- "sim51._AllInts300."
out_abds_all <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

out_abds <- rbind(out_abds, out_abds_all)

CheckAbds(out_abds)
proc_abds <- LabelAbds(out_abds_all)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats)



plCoexist <- PlotCoexistence(pred_stats %>% filter(Mu == -4))
#show(plCoexist)
png("../CavityHOIs-Notes/Coexistence.png", width = 3000, height = 1000, res = 300)
plCoexist
dev.off()

plHist <- PlotHist(proc_abds, c(0.5, 1))
show(plHist)
png("../CavityHOIs-Notes/Histogram.png", width = 3000, height = 1500, res = 300)
plHist
dev.off()

## conceptual understanding plots

mean_concept <- crossing(MuR = 1, Mu = -4, Mean = seq(-0.1, 0.5, length.out = 100),
                         Interaction = c("Pairwise", "Higher Order")) %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order")))

mean_concept <- mean_concept %>% mutate(Solution = ifelse(Interaction == "Pairwise", MuR + Mu * Mean, MuR + Mu * Mean^2))

plMean <- ggplot(mean_concept, aes(x = Mean, y = Solution, color = Interaction)) +
  geom_line(size = 3, alpha = 0.75) + theme_classic() + ylab("Mean") +
  theme(axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15),
        text = element_text(size=20)) + ylim(c(0, 1.1)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.75, linetype = "F1", size = 3)
plMean

tiff("../CavityHOIs-Notes/Concept.tiff", width = 1500, height = 1000, res = 200)
plMean
dev.off()


var_concept <- crossing(Sigma = 0.25, SecAbd = seq(0, 0.2, length.out = 100),
                        Label = c("Identity", "Pairwise", "HOI"))

var_concept <- var_concept %>%
  mutate(Mean = ifelse(Label == "Identity", 0, ifelse(Label == "Pairwise", 0.2, 0.4))) %>%
  mutate(Solution = ifelse(Label == "Identity", SecAbd - Mean^2,
                                                    ifelse(Label == "Pairwise", Mean^2 + Sigma^2 * SecAbd,
                                                           Mean^2 + Sigma^2 * SecAbd^2)))

plVar <- ggplot(var_concept, aes(x = SecAbd, y = Solution, color = Label)) +
  geom_line(size = 1.5, alpha = 0.75) + theme_bw() + ylim(c(0, max(var_concept$SecAbd))) + ylab("SecAbd") +
  theme(legend.title = element_blank()) + ggtitle("Variance Abundance Solution")
plVar


S <- 20
r <- rep(1, times = S)
d <- rep(1, times = S)

mu_A <- -2
sigma_A <- 1
rho_A <- 0
A <- BuildA(S, mu_A, sigma_A, rho_A)

mu_B <- -2
sigma_B <- 1
rho_B <- 0
B <- BuildB(mu_B, sigma_B, rho_B, mu_A, A)

pars <- list(S = S, r = r, d = d, A = A, B = B)

inistate <- runif(S, 0, 0.01)
endtime <- 100
timestep <- 0.1

out_dyn <- IntegrateDynamics(inistate, pars, endtime, timestep, fn = Dynamics)
plSeries <- PlotSeries(out_dyn)
show(plSeries)

tiff("../CavityHOIs-Notes/Dynamics.tiff", width = 1000, height = 1000, res = 300)
plSeries
dev.off()

norm_data <- data.frame(Abundance = seq(-2, 4, length.out = 1000)) %>%
  mutate(Density = dnorm(Abundance, mean = 1, sd = 1))

plDist <- ggplot(norm_data, aes(x = Abundance, y = Density)) +
  geom_line(color = "darkgreen", size = 2) + theme_classic() +
  theme(legend.title = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_blank(),
        strip.text.y = element_text(size = 15),
        text = element_text(size=20)) + ylim(c(0, 0.5)) + xlim(c(0, 4)) +
  xlab("Community Abundance") #+ geom_vline(xintercept = 0, alpha = 0.75, size = 2)
plDist

tiff("../CavityHOIs-Notes/Distribution.tiff", width = 1000, height = 1000, res = 300)
plDist
dev.off()

## Heatmap plotting

GetPredictions <- function(out_data) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredSec = c(), PredFourth = c(),
                     PredV = c(), FuncVal = c(), TermCode = c())
  
  for(row in 1:nrow(out_data)) {
    
    if(row %% 100 == 0) {
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
    ini_var <- ini_sec - ini_mean^2
    
    ini_fourth <- 3 * ini_var^2 + 6 * ini_mean^2 * ini_var + ini_mean^4
    
    h <- function(ini_v, cur_data) with(cur_data, 1 - ini_v * (MuD - 2 * MuB * ini_mean / (S+1) + RhoA * SigmaA^2 * ini_v))
    for(curx in testxs) hxs <- c(hxs, h(curx, cur_data))
    max_v <- testxs[hxs < 0][1]
    ini_v <- uniroot(h, c(0, max_v), cur_data)$root
    
    ini_guess <- c(1 - 0.25 * with(cur_data, sqrt(SigmaA^2 + SigmaB^2)), ini_mean, ini_sec, ini_fourth, ini_v)
    cav_soln <- nleqslv(ini_guess, CavitySoln, method = "Broyden", jac = NULL, cur_data)
    fval <- sum(cav_soln$fvec)
    term_code <- cav_soln$termcd
    cav_soln <- cav_soln$x
    
    max_trials <- 1000
    trial <- 1
    
    while(term_code != 1 && trial < max_trials) {
      ini_guess <- runif(5, 0, 1)
      cav_soln <- nleqslv(ini_guess, CavitySoln, method = "Broyden", jac = NULL, cur_data)
      fval <- sum(cav_soln$fvec)
      term_code <- cav_soln$termcd
      cav_soln <- cav_soln$x
      trial <- trial + 1
    }
    
    PRED <- rbind(PRED, data.frame(PredFraction = ifelse(term_code == 1, cav_soln[1], NA), PredMean = cav_soln[2],
                                   PredSec = cav_soln[3], PredFourth = cav_soln[4], PredV = cav_soln[5],
                                   FuncVal = fval, TermCode = term_code))
  }
  
  ret_out_data <- dplyr::bind_cols(out_data, PRED)
  return(ret_out_data)
}


# some scratch space and code to generate the desired combination of parameters

input_S <- 30
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0#seq(-2, -4, length.out = 50)
input_sigma_A <- 0#seq(0.5, 1, length.out = 2)
input_rho_A <- 0
input_mu_B <- seq(0, -5, length.out = 50)
input_sigma_B <- seq(0.25, 1.5, length.out = 50)
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

#plot_pred <- GetPredictions(LabelAbds(input_params))

plHeatMap <- ggplot(plot_pred, aes(x = Mu, y = Sigma, fill = PredFraction)) +
  geom_tile() + labs(fill = "Phi") +
  scale_fill_viridis(end = 0.9, na.value="darkred") + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme(strip.background = element_rect(fill="white", color = "black"),
                                               aspect.ratio = 1,
                                               panel.background = element_rect(fill = NA),
                                               text = element_text(size=15),
                                               legend.key.size = unit(1, "cm"),
                                               legend.key.width = unit(0.5,"cm")) +
  facet_wrap(~as.factor(Interaction)) + coord_fixed(ratio = (diff(range(plot_pred$Mu)) / diff(range(plot_pred$Sigma)))) +
  labs(x = expression("Mean Interaction Strength"~(mu)), y = expression("Interaction Heterogeneity"~(sigma)),
       fill = expression("Coexisting\nFraction"*(phi)))
show(plHeatMap)
png("../CavityHOIs-Notes/HeatMap.tiff", width = 3000, height = 1000, res = 300)
plHeatMap
dev.off()

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
  fourthN <- x[4]
  v <- x[5]
  
  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * phi^2 * avgN^2))
  var_norm <- with(cur_data, v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + SigmaB^2* phi^2 * secN^2))
  
  
  if(var_norm < 0.001) var_norm <- 0.001
  error_fn <- erf(avg_norm / sqrt(2 * var_norm))
  
  eq1 <- 0.5 * (1 + error_fn)
  eq2 <- avg_norm / 2 * (1 + error_fn) + sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  eq3 <- (avg_norm^2 + var_norm) / 2 * (1 + error_fn)
  eq3 <- eq3 + avg_norm * sqrt(var_norm) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  eq4 <- (3 * var_norm^2 + 6 * avg_norm^2 * var_norm + avg_norm^4) / 2 * (1 + error_fn)
  eq4 <- eq4 + avg_norm * sqrt(var_norm) * (5 * avg_norm * var_norm + avg_norm^3) * exp(- 0.5 * avg_norm^2 / var_norm) / sqrt(2 * pi)
  
  eq1 <- phi - eq1
  eq2 <- avgN - (1 / phi) * eq2
  eq3 <- secN - (1 / phi) * eq3
  eq4 <- fourthN - (1 / phi) * eq4
  eq5 <- with(cur_data, 1 - v * (MuD - 2 * phi * MuB * avgN / S + phi * RhoA * SigmaA^2 * v))
  
  return(c(eq1, eq2, eq3, eq4, eq5))
}

GetStatistics <- function(abds) {
  ret_stats <- abds %>% group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, CommunityID) %>%
    summarise(CommPhi = unique(sum(Abundances > 0) / S), CommMean = mean(Abundances[Abundances > 0]),
              CommSecAbd = mean(Abundances[Abundances > 0]^2), CommVar = CommSecAbd - CommMean^2,
              CommFourthAbd = mean(Abundances[Abundances > 0]^4), Interaction = unique(Interaction)) %>%
    group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB) %>%
    summarise(Phi = mean(CommPhi), MeanAbd = mean(CommMean), SecAbd = mean(CommSecAbd), VarAbd = mean(CommVar),
              FourthAbd = mean(CommFourthAbd), ErrorPhi = sd(CommPhi), ErrorMean = sd(CommMean), ErrorSec = sd(CommSecAbd),
              ErrorVar = sd(CommVar), ErrorFourth = sd(CommFourthAbd), Mu = unique(MuA + MuB),
              Sigma = unique(sqrt(SigmaA^2 + SigmaB^2)), Interaction = unique(Interaction))
  return(ret_stats)
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
pred_stats <- GetPredictions(abd_stats)

pred_stats <- pred_stats %>% mutate(Phi = ifelse(RhoB == 0, Phi, Phi / 2)) %>% mutate(Phi = ifelse(Phi < 0.5, 2 * Phi, Phi))


plCorr <- PlotCorr(pred_stats)
show(plCorr)
tiff("../CavityHOIs-Notes/Correlation.tiff", width = 3000, height = 1000, res = 300)
plCorr
dev.off()


## Multiple Attractor Analysis



