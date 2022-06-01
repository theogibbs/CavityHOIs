
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
  ret_stats <- abds %>% dplyr::group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, CommunityID) %>%
    dplyr::summarise(CommPhi = unique(sum(Abundances > 0) / length(Abundances)), CommMean = mean(Abundances[Abundances > 0]),
                     CommSecAbd = mean(Abundances[Abundances > 0]^2), CommVar = CommSecAbd - CommMean^2,
                     CommFourthAbd = mean(Abundances[Abundances > 0]^4), Interaction = unique(Interaction)) %>%
    dplyr::group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB) %>%
    dplyr::summarise(Phi = mean(CommPhi), MeanAbd = mean(CommMean), SecAbd = mean(CommSecAbd), VarAbd = mean(CommVar),
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

# reading in the data

out_abds <- read.csv("simdata/sim1123_Saturating30.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)


PlotCoexistence <- function(pred_stats) {
  melt_stats <- pred_stats %>%
    ungroup() %>%
    dplyr::select(SigmaR, SigmaA, MuB, SigmaB, Phi, MeanAbd, VarAbd, Interaction) %>%
    melt(id.vars = c("SigmaR", "SigmaA", "MuB", "SigmaB", "Interaction"))
  melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
  melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)
  
  melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                           rep("SAD Mean", times = nrow(pred_stats)),
                           rep("SAD Variance", times = nrow(pred_stats)))
  melt_stats$MuB <- as.factor(melt_stats$MuB)
  melt_stats <- melt_stats %>% filter(variable == "Coexisting Fraction")
  
  plCoexist <- ggplot(melt_stats, aes(x = SigmaB, y = value, color = MuB, shape = MuB)) +
    geom_errorbar(aes(ymin = value - 2 * Error, ymax = value + 2 * Error), alpha = 0.6, width = 0, size = 2) +
    geom_point(size = 7, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                                           text = element_text(size=30),
                                                           strip.text.x = element_text(size = 25),
                                                           strip.text.y = element_text(size = 25),
                                                           legend.text=element_text(size = 25)) +
    geom_line(aes(x = SigmaB, y = Prediction, color = MuB), size = 3, alpha = 0.6) +
    facet_grid(SigmaA~SigmaR, scales = "free",  labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = sigma[A] == .(SigmaA))) +
    labs(x = expression("Interaction Heterogeneity"~(sigma[B])),
         y = expression("Coexisting Fraction"~(phi)),
         color = expression("Interaction\nStrength"(mu[B])),
         shape = expression("Interaction\nStrength"(mu[B])))
  
  return(plCoexist)
}


plHigherOrder <- PlotCoexistence(pred_stats)
plHigherOrder


jpeg("../CavityHOIs-Paper/FigHigherOrder.jpeg", width = 4500, height = 2250, res = 300)
plHigherOrder
dev.off()


