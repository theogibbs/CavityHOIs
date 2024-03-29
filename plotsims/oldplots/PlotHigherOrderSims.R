
source("Functions.R")

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
melt_stats$MuB <- as.factor(melt_stats$MuB)
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
  labs(x = expression("Interaction Heterogeneity"~(sigma[B])),
       y = expression("Coexisting Fraction"~(phi)),
       color = expression("Interaction\nStrength"(mu[B])),
       shape = expression("Interaction\nStrength"(mu[B])))
plSIHigherOrder


jpeg("../CavityHOIs-Paper/SIFigHigherOrder.jpeg", width = 4500, height = 2250, res = 300)
plSIHigherOrder
dev.off()

### now making only cavity method figures

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- c(0, 0.5)
input_rho_A <- 0
input_mu_B <- c(-1, -2, -3)
input_sigma_B <- seq(0, 0.85, length.out = 50)
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_preds <- GetPredictions(input_params, 10)

melt_preds <- out_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, SigmaA, MuB, SigmaB, PredFraction) %>%
  melt(id.vars = c("SigmaR", "SigmaA", "MuB", "SigmaB")) %>%
  mutate(MuB = as.factor(MuB))

plHigherOrder <- ggplot(melt_preds, aes(x = SigmaB, y = value, color = MuB)) +
  geom_line(size = 3, alpha = 0.6) +
  theme_bw() + theme(legend.position = "right",
                     text = element_text(size=30),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 20),
                     legend.text=element_text(size = 25)) +
  facet_grid(SigmaA~SigmaR, scales = "free", labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = sigma[A] == .(SigmaA))) +
  labs(x = expression("Interaction Heterogeneity"~(sigma[B])),
       y = expression("Coexisting Fraction"~(phi)),
       color = expression("Interaction\nStrength"(mu[B]))) +
  ggtitle("(B) Predictions for higher-order interactions")
plHigherOrder

jpeg("../CavityHOIs-Paper/FigHigherOrder.jpeg", width = 4500, height = 2000, res = 300)
plHigherOrder
dev.off()

