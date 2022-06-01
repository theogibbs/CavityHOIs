
source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim1207_Saturating30.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)

melt_stats <- abd_stats %>%
  ungroup() %>%
  dplyr::select(SigmaR, SigmaA, MuB, SigmaB, h, Phi, MeanAbd, VarAbd, Interaction) %>%
  melt(id.vars = c("SigmaR", "SigmaA", "MuB", "SigmaB", "h", "Interaction"))
melt_stats$Error <- c(abd_stats$ErrorPhi, abd_stats$ErrorMean, abd_stats$ErrorVar)

melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(abd_stats)),
                         rep("SAD Mean", times = nrow(abd_stats)),
                         rep("SAD Variance", times = nrow(abd_stats)))
melt_stats$MuB <- as.factor(-melt_stats$MuB)
frac_stats <- melt_stats %>%
  filter(variable == "Coexisting Fraction") %>%
  filter(h == 3)

plSaturating <- ggplot(frac_stats, aes(x = SigmaB, y = value, color = MuB, shape = MuB)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), alpha = 0.6, width = 0, size = 2) +
  geom_point(size = 7, alpha = 0.8) + theme_classic() +
  theme(legend.position = "right",
        panel.spacing = unit(1, "lines"),
        text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        strip.background = element_blank(),
        legend.text=element_text(size = 20)) +
  facet_grid(~SigmaR, scales = "free", labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = h == .(h))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[B])),
       y = expression(atop("Fraction of", "Coexisting Species"~(phi))),
       color = expression("Mean\nInteraction\nStrength"(mu[B])),
       shape = expression("Mean\nInteraction\nStrength"(mu[B])))
plSaturating

jpeg("../CavityHOIs-Paper/figs/FigSaturating.jpeg", width = 3500, height = 1100, res = 300)
plSaturating
dev.off()

ratio_stats <- abd_stats %>%
  mutate(Ratio = MeanAbd / sqrt(VarAbd)) %>%
  ungroup() %>%
  dplyr::select(SigmaR, MuB, SigmaB, h, Ratio) %>%
  melt(id.vars = c("SigmaR", "MuB", "SigmaB", "h")) %>%
  mutate(MuB = as.factor(-MuB))

plSatMeanVar <- ggplot(ratio_stats, aes(x = SigmaB, y = value, color = MuB, shape = MuB)) +
  geom_point(size = 7, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                                         panel.spacing = unit(1, "lines"),
                                                         text = element_text(size=20),
                                                         strip.text.x = element_text(size = 20),
                                                         strip.text.y = element_text(size = 20),
                                                         legend.text=element_text(size = 20)) +
  facet_grid(h~SigmaR, scales = "free",  labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = h == .(h))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[B])),
       y = "Truncated Mean to Standard Deviation Ratio",
       color = expression("Mean\nInteraction\nStrength"(mu[B])),
       shape = expression("Mean\nInteraction\nStrength"(mu[B])))
plSatMeanVar

jpeg("../CavityHOIs-Paper/figs/SIFigSatMeanVar.jpeg", width = 3500, height = 2000, res = 300)
plSatMeanVar
dev.off()



