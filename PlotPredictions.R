#### PACKAGES

library(gridExtra)
library(ggplot2)
library(plyr)
library(tidyverse)
library(nleqslv)
library(viridis)
library(hrbrthemes)
library(cubature)

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


input_S <- 100
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- seq(-2, -5, length.out = 25)
input_sigma_A <- seq(0.25, 2, length.out = 25)
input_rho <- 0
input_mu_B <- 0#c(-2, -3)
input_sigma_B <- 0#seq(0.01, 1, length.out = 50)

abd_cutoff <- 1e-14
gr_cutoff <- 0.001

endtime <- 1e7
inimin <- 0
inimax <- 1

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d, SigmaD = input_sigma_d,
                         MuA = input_mu_A, SigmaA = input_sigma_A, Rho = input_rho, MuB = input_mu_B, SigmaB = input_sigma_B,
                         AbdCutoff = abd_cutoff, GrCutoff = gr_cutoff, EndTime = endtime, IniMin = inimin, IniMax = inimax)

temp_params <- input_params

input_mu_B <- input_mu_A
input_sigma_B <- input_sigma_A
input_mu_A <- 0
input_sigma_A <- 0

input_params <- rbind(input_params,
                      crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d, SigmaD = input_sigma_d,
                         MuA = input_mu_A, SigmaA = input_sigma_A, Rho = input_rho, MuB = input_mu_B, SigmaB = input_sigma_B,
                         AbdCutoff = abd_cutoff, GrCutoff = gr_cutoff, EndTime = endtime, IniMin = inimin, IniMax = inimax))

temp_params$MuA <- temp_params$MuA / 2
temp_params$MuB<- temp_params$MuA
temp_params$SigmaA <- temp_params$SigmaA / sqrt(2)
temp_params$SigmaB <- temp_params$SigmaA

input_params <- rbind(input_params, temp_params)

plot_pred <- get_preds(input_params)

plFracA <- ggplot(plot_pred %>% filter(SigmaB == 0), aes(x = SigmaA, y = PredFraction, color = MuA)) + 
  geom_point() + theme_bw() + labs(color = "MuA") + #xlim(c(0.15, 0.75)) + ylim(c(0.825, 1)) +
  ggtitle("Pairwise")

plFracB <- ggplot(plot_pred %>% filter(SigmaB != 0), aes(x = Sigma, y = PredFraction, color = Mu)) + 
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
  geom_raster(interpolate = TRUE) + labs(fill = "Phi") + scale_fill_viridis(end = 0.9) + scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + theme(strip.background = element_rect(fill="white", color = "black"),
                                               aspect.ratio = 1,
                                               panel.background = element_rect(fill = NA),
                                               text = element_text(size=15),
                                               legend.key.size = unit(1, "cm"),
                                               legend.key.width = unit(0.5,"cm")) +
  facet_wrap(~as.factor(PlotType)) + coord_fixed(ratio = (diff(range(plot_pred$Mu)) / diff(range(plot_pred$Sigma))))
plHeatMap

