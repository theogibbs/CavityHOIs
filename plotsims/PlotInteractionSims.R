
source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim1202_Pairwise30.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)

melt_stats <- pred_stats %>%
  ungroup() %>%
  dplyr::select(SigmaR, Mu, Sigma, Phi, MeanAbd, VarAbd, Interaction) %>%
  melt(id.vars = c("SigmaR", "Mu", "Sigma", "Interaction"))
melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)

melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                         rep("SAD Mean", times = nrow(pred_stats)),
                         rep("SAD Variance", times = nrow(pred_stats)))
melt_stats$Mu <- as.factor(-melt_stats$Mu)

plSIPairwise <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 2) +
  geom_point(size = 7) + theme_bw() + theme(legend.position = "right",
                                            text = element_text(size=30),
                                            strip.text.x = element_text(size = 25),
                                            strip.text.y = element_text(size = 20),
                                            legend.text=element_text(size = 25)) +
  geom_line(aes(x = Sigma, y = Prediction, color = Mu), size = 3, alpha = 0.6) +
  facet_grid(variable~SigmaR, switch = "y", scales = "free", labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = .(variable))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A])),
       y = "",
       color = expression("Interaction\nStrength"(mu[A])),
       shape = expression("Interaction\nStrength"(mu[A]))) +
  scale_y_continuous(position = "right")
plSIPairwise

jpeg("../CavityHOIs-Paper/figs/SIFigPairwise.jpeg", width = 4500, height = 3000, res = 300)
plSIPairwise
dev.off()

# reading in the data

out_abds <- read.csv("simdata/sim1125_HigherOrder30.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)

melt_stats <- pred_stats %>%
  ungroup() %>%
  dplyr::select(SigmaR, SigmaA, MuB, SigmaB, Phi, MeanAbd, VarAbd, Interaction) %>%
  melt(id.vars = c("SigmaR", "SigmaA", "MuB", "SigmaB", "Interaction"))
melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
melt_stats$Prediction <- c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)

melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                         rep("SAD Mean", times = nrow(pred_stats)),
                         rep("SAD Variance", times = nrow(pred_stats)))
melt_stats$MuB <- as.factor(- melt_stats$MuB)
melt_stats <- melt_stats %>% filter(variable == "Coexisting Fraction")

plSIHigherOrder <- ggplot(melt_stats, aes(x = SigmaB, y = value, color = MuB, shape = MuB)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), alpha = 0.6, width = 0, size = 2) +
  geom_point(size = 7, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                                         text = element_text(size=30),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size = 25)) +
  geom_line(aes(x = SigmaB, y = Prediction, color = MuB), size = 3, alpha = 0.6) +
  facet_grid(SigmaA~SigmaR, scales = "free",  labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = sigma[A] == .(SigmaA))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[B])),
       y = expression("Coexisting Fraction"~(phi)),
       color = expression("Interaction\nStrength"(mu[B])),
       shape = expression("Interaction\nStrength"(mu[B])))
plSIHigherOrder

jpeg("../CavityHOIs-Paper/figs/SIFigHigherOrder.jpeg", width = 4500, height = 2250, res = 300)
plSIHigherOrder
dev.off()

### now making only cavity method figures

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- c(-1, -2, -3)
input_sigma_A <- seq(0.01, 0.85, length.out = 50)
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_pw_preds <- GetPredictions(input_params, 10)

melt_pw_preds <- out_pw_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuA, SigmaA, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuA", "SigmaA")) %>%
  mutate(MuA = as.factor(-MuA)) %>%
  mutate(Interaction = "(A) Pairwise") %>%
  dplyr::rename(Sigma = SigmaA, Mu = MuA)

### higher-order interactions

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0
input_rho_A <- 0
input_mu_B <- c(-1, -2, -3)
input_sigma_B <- seq(0, 0.85, length.out = 50)
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_hoi_preds <- GetPredictions(input_params, 10)

melt_hoi_preds <- out_hoi_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuB, SigmaB, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuB", "SigmaB")) %>%
  mutate(MuB = as.factor(-MuB)) %>%
  mutate(Interaction = "(B) Higher-order") %>%
  dplyr::rename(Sigma = SigmaB, Mu = MuB)

# mixed interactions

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- 0.5
input_rho_A <- 0
input_mu_B <- c(-1, -2, -3)
input_sigma_B <- seq(0, 0.85, length.out = 50)
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_mx_preds <- GetPredictions(input_params, 10)

melt_mx_preds <- out_mx_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuB, SigmaB, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuB", "SigmaB")) %>%
  mutate(MuB = as.factor(-MuB)) %>%
  mutate(Interaction = "(C) Pairwise and higher-order") %>%
  dplyr::rename(Sigma = SigmaB, Mu = MuB)

melt_preds <- rbind(melt_pw_preds, melt_hoi_preds, melt_mx_preds)

plInteractions <- ggplot(melt_preds, aes(x = Sigma, y = value, group = Mu)) +
  geom_line(size = 3, alpha = 0.6, aes(color = Mu, linetype = Mu)) +
  theme_classic() + theme(legend.position = "top",
                     text = element_text(size=30),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 25),
                     strip.background = element_blank(),
                     legend.text=element_text(size = 25)) +
  facet_grid(Interaction ~ SigmaR,
             labeller = label_bquote(rows = .(Interaction), cols = sigma[R] == .(SigmaR)),
             switch = "y", scales = "free") +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A]~"or"~sigma[B])),
       y = expression("Fraction of Coexisting Species"~(phi)),
       color = expression("Mean Interaction Strength"~(mu[A]~"or"~mu[B])),
       linetype = expression("Mean Interaction Strength"~(mu[A]~"or"~mu[B]))) +
  scale_y_continuous(position = "right") +
  scale_linetype_manual(values=c(1, 2, 6)) +
  guides(linetype=guide_legend(keywidth = 8, keyheight = 1),
         color=guide_legend(keywidth = 8, keyheight = 1))
plInteractions

jpeg("../CavityHOIs-Paper/figs/FigInteractions.jpeg", width = 6000, height = 5000, res = 300)
plInteractions
dev.off()

# mean interaction strength figure

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- seq(-2, -6, length.out = 50)
input_sigma_A <- c(0.1, 0.5, 1)
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_pw_preds <- GetPredictions(input_params, 10)

melt_pw_preds <- out_pw_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuA, SigmaA, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuA", "SigmaA")) %>%
  mutate(MuA = -MuA) %>%
  mutate(Interaction = "(A) Pairwise") %>%
  dplyr::rename(Sigma = SigmaA, Mu = MuA)

### higher-order interactions

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0
input_sigma_A <- 0
input_rho_A <- 0
input_mu_B <- seq(-2, -6, length.out = 50)
input_sigma_B <- c(0.1, 0.5, 1)
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_hoi_preds <- GetPredictions(input_params, 10)

melt_hoi_preds <- out_hoi_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuB, SigmaB, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuB", "SigmaB")) %>%
  mutate(MuB = -MuB) %>%
  mutate(Interaction = "(B) Higher-order") %>%
  dplyr::rename(Sigma = SigmaB, Mu = MuB)

# mixed interactions

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- 0.5
input_rho_A <- 0
input_mu_B <- seq(-2, -6, length.out = 50)
input_sigma_B <- c(0.1, 0.5, 1)
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_mx_preds <- GetPredictions(input_params, 10)

melt_mx_preds <- out_mx_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuB, SigmaB, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuB", "SigmaB")) %>%
  mutate(MuB = -MuB) %>%
  mutate(Interaction = "(C) Pairwise and higher-order") %>%
  dplyr::rename(Sigma = SigmaB, Mu = MuB)

melt_preds <- rbind(melt_pw_preds, melt_hoi_preds, melt_mx_preds) %>%
  mutate(Sigma = as.factor(Sigma))

plMeanInteractions <- ggplot(melt_preds, aes(x = Mu, y = value, group = Sigma)) +
  geom_line(size = 3, alpha = 0.6, aes(color = Sigma, linetype = Sigma)) +
  theme_classic() + theme(legend.position = "top",
                          text = element_text(size=30),
                          strip.text.x = element_text(size = 25),
                          strip.text.y = element_text(size = 25),
                          strip.background = element_blank(),
                          legend.text=element_text(size = 25)) +
  facet_grid(Interaction ~ SigmaR,
             labeller = label_bquote(rows = .(Interaction), cols = sigma[R] == .(SigmaR)),
             switch = "y", scales = "free") +
  labs(x = expression("Mean Interaction Strength"~(mu[A]~"or"~mu[B])),
       y = expression("Fraction of Coexisting Species"~(phi)),
       color = expression("Variation in Interaction Strengths"~(sigma[A]~"or"~sigma[B])),
       linetype = expression("Variation in Interaction Strengths"~(sigma[A]~"or"~sigma[B]))) +
  scale_y_continuous(position = "right") +
    scale_linetype_manual(values=c(1, 2, 6)) +
    guides(linetype=guide_legend(keywidth = 8, keyheight = 1),
           color=guide_legend(keywidth = 8, keyheight = 1))
plMeanInteractions

jpeg("../CavityHOIs-Paper/figs/SIFigMeanInts.jpeg", width = 6000, height = 5000, res = 300)
plMeanInteractions
dev.off()
