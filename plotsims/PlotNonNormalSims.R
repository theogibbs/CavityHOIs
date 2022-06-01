
source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim_NonNormal.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)

melt_stats <- pred_stats %>%
  ungroup() %>%
  dplyr::select(S, Mu, Sigma, Phi, Interaction) %>%
  melt(id.vars = c("S", "Mu", "Sigma", "Interaction")) %>%
  mutate(Error = pred_stats$ErrorPhi) %>%
  mutate(Prediction = pred_stats$PredFraction) %>%
  mutate(Mu = as.factor(-Mu))

plNonNormal <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 2) +
  geom_point(size = 7) + theme_classic() +
  theme(legend.position = "right",
        panel.spacing = unit(1.5, "lines"),
        text = element_text(size=30),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        strip.background = element_blank(),
        legend.text=element_text(size = 25)) +
  geom_line(aes(x = Sigma, y = Prediction, color = Mu), size = 3, alpha = 0.6) +
  facet_grid(~S, scales = "free", labeller = label_bquote(cols = .(S) ~ "Species")) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[B])),
       y = expression(atop("Fraction of","Coexisting Species"~(phi))),
       color = expression(atop("Mean\nInteraction\nStrength", (mu[B]))),
       shape = expression(atop("Mean\nInteraction\nStrength", (mu[B])))) +
  ggtitle("Uniformly distributed growth rates and interactions")
plNonNormal

jpeg("../CavityHOIs-Paper/figs/SIFigNonNormal.jpeg", width = 4500, height = 1500, res = 300)
plNonNormal
dev.off()

