
### now making only cavity method figures

input_S <- 20
input_mu_r <- 1.5
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_sigma_A <- seq(0.001, 0.5, length.out = 100) * sqrt(input_S)
input_mu_A <- 0
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)
input_params$MuA <- - sqrt(3 * input_S) * input_sigma_A

input_sigma_B <- input_sigma_A * sqrt(input_S)
input_mu_B <- 0
input_mu_A <- 0
input_sigma_A <- 0

new_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                       SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                       MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h)

new_params$MuB <- - sqrt(3) * input_S * input_sigma_B

input_params <- rbind(input_params, new_params)

out_preds <- GetPredictions(input_params, 10)
out_preds <- LabelAbds(out_preds)

melt_preds <- out_preds %>%
  ungroup() %>%
  mutate(PredSD = PredSec - PredMean^2) %>%
  filter(PredSD < 1) %>%
  dplyr::select(Mu, Sigma, PredFraction, PredMean, PredSD, Interaction) %>%
  melt(id.vars = c("Mu", "Sigma", "Interaction")) %>%
  filter(variable == "PredFraction")

plConstrained <- ggplot(melt_preds, aes(x = - 2 * Mu, y = value)) +
  geom_line(size = 3, alpha = 0.6) +
  theme_bw() + theme(legend.position = "right",
                     text = element_text(size=30),
                     strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 20),
                     legend.text=element_text(size = 25)) +
  facet_wrap(~Interaction, scales = "free") +
  labs(x = "Maximum Interaction Strength",
       y = expression(atop("Fraction of", "Coexisting Species"~(phi))))
plConstrained

jpeg("../CavityHOIs-Paper/figs/SIFigConstrained.jpeg", width = 3500, height = 1500, res = 300)
plConstrained
dev.off()
