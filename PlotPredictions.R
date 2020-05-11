#### PACKAGES

library(gridExtra)
library(nleqslv)
library(VGAM)

source("Functions.R")

cavity_soln <- function(x, cur_data) {
  
  phi <- x[1]
  avgN <- x[2]
  secN <- x[3]
  fourthN <- x[4]
  v <- x[5]
  
  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * (phi^2 * avgN^2 * (S-2) / (S+1) + phi * secN / (S+1))))
  var_norm <- abs(with(cur_data,
                       v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + SigmaB^2 * ((S-2) * phi * secN^2 + fourthN) * phi / (S+1))))
  
  f2 <- function(N, avg_norm, var_norm) {
    return(N * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }
  f3 <- function(N, avg_norm, var_norm) {
    return(N^2 * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }
  f4 <- function(N, avg_norm, var_norm) {
    return(N^4 * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }
  
  eq1 <- phi - (1 - pnorm(0, mean = avg_norm, sd = sqrt(var_norm)))
  eq2 <- avgN - (1 / phi) * cubintegrate(f2, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq3 <- secN - (1 / phi) * cubintegrate(f3, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq4 <- fourthN - (1 / phi) * cubintegrate(f4, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq5 <- with(cur_data, 1 - v * (MuD - 2 * phi * MuB * avgN / (S+1) + phi * Rho * SigmaA^2 * v))
  
  return(c(eq1, eq2, eq3, eq4, eq5))
}
get_preds <- function(out_data) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredVar = c())
  
  for(row in 1:nrow(out_data)) {
    print("Row number:"); print(row)
    
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
    
    h <- function(ini_v, cur_data) with(cur_data, 1 - ini_v * (MuD - 2 * MuB * ini_mean / (S+1) + Rho * SigmaA^2 * ini_v))
    for(curx in testxs) hxs <- c(hxs, h(curx, cur_data))
    max_v <- testxs[hxs < 0][1]
    ini_v <- uniroot(h, c(0, max_v), cur_data)$root
    
    ini_guess <- c(1, ini_mean, ini_sec, ini_fourth, ini_v)
    print(ini_guess)
    
    cav_soln <- nleqslv(ini_guess, cavity_soln, method = "Broyden", jac = NULL, cur_data)$x
    print(cav_soln)
    PRED <- rbind(PRED, data.frame(PredFraction = cav_soln[1], PredMean = cav_soln[2],
                                   PredVar = cav_soln[3], PredFourth = cav_soln[4]))
  }
  ret_out_data <- dplyr::bind_cols(out_data, PRED)
  return(ret_out_data)
}



####

filename <- "sim501_OnlyHOIs300.csv"
filename <- paste0("simdata/", filename)
out_data <- read.csv(filename, header = TRUE, sep = ",")
out_data <- rbind(out_data, read.csv("simdata/sim503_OnlyHOIs300.csv", header = TRUE, sep = ","))
out_data <- rbind(out_data, read.csv("simdata/sim504_OnlyHOIs300.csv", header = TRUE, sep = ","))
#out_data <- rbind(out_data, read.csv("simdata/sim505_OnlyHOIs300.csv", header = TRUE, sep = ","))


out_data %>% filter(Equilibrium == FALSE) %>%
  summarise(NonEqNumRuns = length(Equilibrium))

out_data %>% filter(Equilibrium == TRUE, Uninvasible == FALSE) %>%
  summarise(InvNumRuns = length(Equilibrium))

plot_coexist <- out_data %>% filter(Equilibrium == TRUE)
plot_coexist$FracCoexist <- plot_coexist$NumCoexist / plot_coexist$S

plot_coexist <- summarySE(plot_coexist, measurevar = "FracCoexist",
                          groupvars = c("S", "MuR", "SigmaR", "MuD", "SigmaD", "MuA", "SigmaA",
                                        "Rho", "MuB", "SigmaB"))
plot_coexist <- get_preds(plot_coexist)

plCoexist <- ggplot(plot_coexist, aes(x = SigmaB, y = FracCoexist, color = as.factor(S))) + 
  geom_errorbar(aes(ymin= FracCoexist - se, ymax = FracCoexist + se), width = 0) +
  geom_line(aes(x = SigmaB, y = PredFraction, color = as.factor(S)), size = 1.5, alpha = 0.75) +
  geom_point(size = 2.5) + theme_bw() + labs(color = "S") + #xlim(c(0.15, 0.75)) + ylim(c(0.825, 1)) +
  ggtitle("Coexistence Predictions")
plCoexist


plot_mean <- out_data
plot_mean <- summarySE(plot_mean, measurevar = "MeanAbd",
                          groupvars = c("S", "MuR", "SigmaR", "MuD", "SigmaD", "MuA", "SigmaA",
                                        "Rho", "MuB", "SigmaB"))
plot_mean <- get_preds(plot_mean)

plMean <- ggplot(plot_mean, aes(x = SigmaA, y = MeanAbd, color = MuB)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin= MeanAbd - se, ymax = MeanAbd + se), width = 0.02) +
  labs(color = "MuB") + ylab("Mean Abundance") +
  geom_line(aes(x = SigmaA, y = PredMean), size = 2.5, alpha = 0.5) +
  theme_bw() + ggtitle("Mean Abundance")
plMean


### only predictions

cavity_soln <- function(x, cur_data) {
  
  phi <- x[1]
  avgN <- x[2]
  secN <- x[3]
  fourthN <- x[4]
  v <- x[5]

  avg_norm <- with(cur_data, v * (MuR + phi * MuA * avgN +  MuB * (phi^2 * avgN^2 * (S-2) / (S+1) + phi * secN / (S+1))))
  var_norm <- abs(with(cur_data,
                       v^2 * (SigmaR^2 + SigmaA^2 * phi * secN + SigmaB^2 * ((S-2) * phi * secN^2 + fourthN) * phi / (S+1))))

    f2 <- function(N, avg_norm, var_norm) {
    return(N * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }
  f3 <- function(N, avg_norm, var_norm) {
    return(N^2 * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }
  f4 <- function(N, avg_norm, var_norm) {
    return(N^4 * dnorm(N, mean = avg_norm, sd = sqrt(var_norm)))
  }

  eq1 <- phi - (1 - pnorm(0, mean = avg_norm, sd = sqrt(var_norm)))
  eq2 <- avgN - (1 / phi) * cubintegrate(f2, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq3 <- secN - (1 / phi) * cubintegrate(f3, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq4 <- fourthN - (1 / phi) * cubintegrate(f4, lower = 0, upper = Inf, avg_norm = avg_norm, var_norm = var_norm)$integral
  eq5 <- with(cur_data, 1 - v * (MuD - 2 * phi * MuB * avgN / (S+1) + phi * Rho * SigmaA^2 * v))
  
  return(c(eq1, eq2, eq3, eq4, eq5))
}

get_preds <- function(out_data) {
  PRED <- data.frame(PredFraction = c(), PredMean = c(), PredVar = c(), PredFourth = c(), PredV = c())
  
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
    
    h <- function(ini_v, cur_data) with(cur_data, 1 - ini_v * (MuD - 2 * MuB * ini_mean / (S+1) + Rho * SigmaA^2 * ini_v))
    for(curx in testxs) hxs <- c(hxs, h(curx, cur_data))
    max_v <- testxs[hxs < 0][1]
    ini_v <- uniroot(h, c(0, max_v), cur_data)$root
    
    ini_guess <- c(1 - 0.25 * with(cur_data, sqrt(SigmaA^2 + SigmaB^2)), ini_mean, ini_sec, ini_fourth, ini_v)
    ini_guess[1] <- 0.1
    
    new_cav_soln <- nleqslv(ini_guess, cavity_soln, method = "Broyden", jac = NULL, cur_data)$x
    ini_guess[1] <- 0.1
    
    if((new_cav_soln[1] < cur_data$GrCutoff) || (new_cav_soln[1] > 1.0001)) {
      ini_guess[1] <- 0.1
      new_cav_soln <- nleqslv(ini_guess, cavity_soln, method = "Broyden", jac = NULL, cur_data)$x
      
    }
    
    cav_soln <- new_cav_soln
    
    cur_phi <- cav_soln[1]
    cur_sec <- cav_soln[3]
    cur_v <- cav_soln[5]
    
    #if(with(cur_data, cur_phi * SigmaA^2 + cur_phi^2 * SigmaB^2 > 1)) {
    #  cav_soln <- rep(0, times = 5)
    #}
    
    #print(test_thresh)
    #if(test_thresh < 0) {
    #  cav_soln <- rep(0, times = 5)
    #}
    
    PRED <- rbind(PRED, data.frame(PredFraction = cav_soln[1], PredMean = cav_soln[2],
                                   PredVar = cav_soln[3], PredFourth = cav_soln[4], PredV = cav_soln[5]))
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
input_sigma_B <- seq(0.5, 1.5, length.out = 50)
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

plot_pred <- GetPredictions(input_params)

plFracA <- ggplot(plot_pred %>% filter(SigmaB == 0), aes(x = SigmaA, y = PredFraction, color = MuA)) + 
  geom_point() + theme_bw() + labs(color = "MuA") + #xlim(c(0.15, 0.75)) + ylim(c(0.825, 1)) +
  ggtitle("Pairwise")

plFracB <- ggplot(plot_pred %>% filter(SigmaB != 0), aes(x = SigmaB, y = PredFraction, color = MuB)) + 
  geom_point() + theme_bw() + labs(color = "Mu") + #xlim(c(0.15, 0.75)) + ylim(c(0.825, 1)) +
  ggtitle("HOI")
grid.arrange(plFracA, plFracB, nrow = 1)

plot_pred <- plot_pred %>%
  mutate(InteractionType = ifelse(SigmaA != 0, ifelse(SigmaB == 0, "Pairwise Interactions", "Mixed Interactions"),
                                  "Higher Order Interactions")) %>%
  mutate(PlotType = factor(InteractionType,
                           levels = c("Pairwise Interactions", "Mixed Interactions", "Higher Order Interactions"))) %>%
  mutate(Sigma = sqrt(SigmaA^2 + SigmaB^2)) %>%
  mutate(Mu = MuA + MuB)

plHeatMap <- ggplot(plot_pred, aes(x = Mu, y = Sigma, fill = PredFraction)) +
  geom_raster(interpolate = TRUE) + labs(fill = "Phi") +
  scale_fill_viridis(end = 0.9) + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme(strip.background = element_rect(fill="white", color = "black"),
                                               aspect.ratio = 1,
                                               panel.background = element_rect(fill = NA),
                                               text = element_text(size=15),
                                               legend.key.size = unit(1, "cm"),
                                               legend.key.width = unit(0.5,"cm")) +
  facet_wrap(~as.factor(PlotType)) + coord_fixed(ratio = (diff(range(plot_pred$Mu)) / diff(range(plot_pred$Sigma))))
plHeatMap




###### New plotting code

### only predictions

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
  eq5 <- with(cur_data, 1 - v * (MuD - 2 * phi * MuB * avgN / (S+1) + phi * RhoA * SigmaA^2 * v))
  
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

LabelAbds <- function(abds) {
  
  EqRuns <- abds %>%
    filter(Equilibrium == FALSE) %>%
    summarise(NonEqNumRuns = length(Equilibrium))
  
  if(EqRuns$NonEqNumRuns != 0) print("WARNING: Some runs did not reach equilibrium.")
  
  InvRuns <- abds %>%
    filter(Equilibrium == TRUE, Uninvadeable == FALSE) %>%
    summarise(InvNumRuns = length(Equilibrium))
  
  if(InvRuns$InvNumRuns != 0) print("WARNING: Some equilibria are invasible.")
  
  ret_abds <- abds %>%
    mutate(Interaction= ifelse(SigmaA != 0, ifelse(SigmaB == 0, "Pairwise", "Mixed"), "Higher Order")) %>%
    mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Mixed", "Higher Order"))) %>%
    mutate(Sigma = signif(sqrt(SigmaA^2 + SigmaB^2)))
  
  return(ret_abds)
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

PlotCoexistence <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    select(Mu, Sigma, Phi, MeanAbd, VarAbd, Interaction) %>%
    melt(id.vars = c("Mu", "Sigma", "Interaction"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)), rep("SAD Mean", times = nrow(pred_stats)),
                           rep("SAD Variance", times = nrow(pred_stats)))
  
  plCoexist <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Interaction, shape = Interaction)) +
    geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 1) +
    geom_point(size = 4) + theme_bw() + theme(axis.title.y=element_blank(),
                                              strip.text.x = element_text(size = 15),
                                              strip.text.y = element_text(size = 15),
                                              text = element_text(size=20)) +
    geom_line(aes(x = Sigma, y = Prediction, color = Interaction), size = 2, alpha = 0.75) +
    facet_wrap(~ variable, scales = "free")
    
  return(plCoexist)
}

PlotHist <- function(proc_abds, hist_sigmas) {
  plot_abds <- proc_abds %>%
    filter(Sigma == hist_sigmas[1] | Sigma == hist_sigmas[2]) %>%
    filter(Abundances > 0)
  
  pred_data <- plot_abds %>%
    select(-one_of("Abundances", "SpeciesID", "CommunityID", "Index")) %>%
    unique()
  
  pred_abds <- GetPredictions(pred_data)
  plot_abds <- merge(plot_abds, pred_abds)
  plot_abds <- plot_abds %>%
    mutate(PredDensity = dnorm(Abundances, mean = PredMean, sd = sqrt(PredSec - PredMean^2)))
  
  plHist <- ggplot(plot_abds, aes(x = Abundances, y = ..density..)) +
    geom_histogram(fill = "white", color = "black", bins = 40) +
    ggtitle("Species Abundance Distributions") +
    facet_grid(Sigma ~ Interaction, labeller = label_bquote(rows = sigma == .(Sigma))) +
    theme_bw() + theme(strip.text.x = element_text(size = 15),
                       strip.text.y = element_text(size = 15),
                       text = element_text(size=20)) + ylab("Density") +
    geom_line(aes(x = Abundances, y = PredDensity), color = "blue", size = 1)
  return(plHist)
}

fileregex <- "sim508_OnlyPairwise300."
out_abds_pw <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))
fileregex <- "sim510_OnlyHOIs300."
out_abds_hois <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))
out_abds <- rbind(out_abds_pw, out_abds_hois)


proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats)


plCoexist <- PlotCoexistence(pred_stats)
show(plCoexist)
tiff("../CavityHOIs-Notes/Coexistence.tiff", width = 3000, height = 1000, res = 200)
plCoexist
dev.off()

plHist <- PlotHist(proc_abds, c(0.5, 1))
show(plHist)
tiff("../CavityHOIs-Notes/Histogram.tiff", width = 3000, height = 1000, res = 200)
plHist
dev.off()

## conceptual understanding

conc_data <- crossing(MuR = 1, Mu = -4, Mean = seq(-0.1, 0.5, length.out = 100), Label = c("Identity", "Pairwise", "HOI"))

conc_data <- conc_data %>% mutate(Solution = ifelse(Label == "Identity", Mean,
                                                    ifelse(Label == "Pairwise", MuR + Mu * Mean, MuR + Mu * Mean^2)))

ggplot(conc_data, aes(x = Mean, y = Solution, color = Label)) +
  geom_line(size = 1.5, alpha = 0.75) + theme_bw() + ylim(c(0, 1)) + ylab("Mean") +
  theme(legend.title = element_blank()) + ggtitle("Mean Abundance Solution")



