
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
melt_stats$Mu <- as.factor(melt_stats$Mu)

plSIPairwise <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 2) +
  geom_point(size = 7) + theme_bw() + theme(legend.position = "right",
                                            text = element_text(size=30),
                                            strip.text.x = element_text(size = 25),
                                            strip.text.y = element_text(size = 20),
                                            legend.text=element_text(size = 25)) +
  geom_line(aes(x = Sigma, y = Prediction, color = Mu), size = 3, alpha = 0.6) +
  facet_grid(variable~SigmaR, switch = "y", scales = "free", labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = .(variable))) +
  labs(x = expression("Interaction Heterogeneity"~(sigma[A])),
       y = "",
       color = expression("Interaction\nStrength"(mu[A])),
       shape = expression("Interaction\nStrength"(mu[A]))) +
  scale_y_continuous(position = "right")
plSIPairwise

jpeg("../CavityHOIs-Paper/SIFigPairwise.jpeg", width = 4500, height = 3000, res = 300)
plSIPairwise
dev.off()

### now making only cavity method figures

input_S <- 30
input_mu_r <- 1.5
input_sigma_r <- c(0, 0.25, 0.5)
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- c(0, -1, -2)
input_sigma_A <- seq(0.01, 0.85, length.out = 50)
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

out_preds <- GetPredictions(input_params, 10)

melt_preds <- out_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuA, SigmaA, PredFraction) %>%
  melt(id.vars = c("SigmaR", "MuA", "SigmaA")) %>%
  mutate(MuA = as.factor(MuA))

plPairwise <- ggplot(melt_preds, aes(x = SigmaA, y = value, color = MuA)) +
  geom_line(size = 3, alpha = 0.6) +
  theme_bw() + theme(legend.position = "right",
                     text = element_text(size=30),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 20),
                     legend.text=element_text(size = 25)) +
  facet_wrap(~SigmaR, labeller = label_bquote(cols = sigma[R] == .(SigmaR))) +
  labs(x = expression("Interaction Heterogeneity"~(sigma[A])),
       y = expression("Coexisting Fraction"~(phi)),
       color = expression("Interaction\nStrength"(mu[A])),
       shape = expression("Interaction\nStrength"(mu[A]))) +
  ggtitle("(A) Predictions for pairwise interactions")
plPairwise

jpeg("../CavityHOIs-Paper/FigPairwise.jpeg", width = 4500, height = 1500, res = 300)
plPairwise
dev.off()
