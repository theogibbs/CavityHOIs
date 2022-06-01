library(viridis)

# some scratch space and code to generate the desired combination of parameters

input_S <- 300
input_mu_r <- 1
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- 0#seq(-2, -4, length.out = 50)
input_sigma_A <- 0#seq(0.5, 1, length.out = 2)
input_rho_A <- c(-0.5, 0, 0.5)
input_mu_B <- seq(0, -5, length.out = 50)
input_sigma_B <- seq(0.01, 1.5, length.out = 50)
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

plHeatMap <- ggplot(plot_pred, aes(x = -Mu, y = Sigma, fill = PredFraction)) +
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
  labs(x = expression("Mean Interaction Strength"~(mu[A] + mu[B])),
       y = expression("Variability in Interaction Strengths"~(sqrt(sigma[A]^2 + sigma[B]^2))),
       fill = expression("Fraction of\nCoexisting\nSpecies"*(phi)))
show(plHeatMap)

jpeg("../CavityHOIs-Paper/figs/SIFigHeatMap.jpeg", width = 3000, height = 3000, res = 300)
plHeatMap
dev.off()
