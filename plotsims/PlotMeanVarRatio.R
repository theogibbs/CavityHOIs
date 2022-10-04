
source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim1125_Pairwise30.csv", header = TRUE)
out_abds <- rbind(out_abds, read.csv("simdata/sim1125_HigherOrder30.csv", header = TRUE))

CheckAbds(out_abds)
out_abds <- out_abds %>%
  filter(Equilibrium == TRUE, Uninvadeable == TRUE)

proc_abds <- LabelAbds(out_abds) %>%
  filter(Interaction != "Mixed") %>%
  filter(SigmaR != 0.5) %>%
  filter(Sigma > 0.01) %>%
  filter(!(Mu == -3 & Interaction == "Pairwise")) %>%
  filter(!(Mu == -4 & Interaction == "Pairwise")) %>%
  filter(SigmaA != 0.5)

abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)
adj_preds <- GetAdjStats(pred_stats) %>%
  mutate(Ratio = AdjMean / AdjSD) %>%
  mutate(PredRatio = AdjPredMean / AdjPredSD)

melt_stats <- adj_preds %>%
  ungroup() %>%
  dplyr::select(SigmaR, Mu, Sigma, Ratio, Interaction) %>%
  melt(id.vars = c("SigmaR", "Mu", "Sigma", "Interaction")) %>%
  mutate(Prediction = adj_preds$PredRatio) %>%
  mutate(Mu = as.factor(-Mu)) %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order")))

plSIMeanVar <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
  geom_point(size = 7, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                                         text = element_text(size=30),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size = 25)) +
  geom_line(aes(x = Sigma, y = Prediction, color = Mu), size = 1, alpha = 0.6) +
  facet_grid(SigmaR~Interaction, scales = "free",  labeller = label_bquote(cols = .(Interaction), rows = sigma[R] == .(SigmaR))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A]  ~ "or" ~ sigma[B])),
       y = expression(atop("Non-truncated Mean to", "Standard Deviation Ratio"~(mu[0]/sigma[0]))),
       color = expression("Mean\nInteraction\nStrength"(mu[A]  ~ "or" ~ mu[B])),
       shape = expression("Mean\nInteraction\nStrength"(mu[A]  ~ "or" ~ mu[B])))
plSIMeanVar

jpeg("../CavityHOIs-Paper/figs/SIFigMeanVar.jpeg", width = 4500, height = 3000, res = 300)
plSIMeanVar
dev.off()

