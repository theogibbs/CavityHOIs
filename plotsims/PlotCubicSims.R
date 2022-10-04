## Requirements

source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim0114_CubicHOIs.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds) %>%
  mutate(SelfReg = ifelse(MuD == 1, "Constant", "Mean-Adjusted"))
pred_stats <- abd_stats

melt_stats <- pred_stats %>%
  ungroup() %>%
  dplyr::select(SigmaR, SigmaA, MuB, SigmaB, SelfReg, Phi, MeanAbd, VarAbd, Interaction) %>%
  melt(id.vars = c("SigmaR", "SigmaA", "MuB", "SigmaB", "SelfReg", "Interaction"))
melt_stats$Error <- c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)
melt_stats$variable <- c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                         rep("SAD Mean", times = nrow(pred_stats)),
                         rep("SAD Variance", times = nrow(pred_stats)))
melt_stats$MuB <- as.factor(- melt_stats$MuB)
melt_stats <- melt_stats %>% filter(variable == "Coexisting Fraction")

plSICubic <- ggplot(melt_stats, aes(x = SigmaB, y = value, color = MuB, shape = MuB)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), alpha = 0.6, width = 0, size = 2) +
  geom_point(size = 7, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                                         panel.spacing = unit(1.5, "lines"),
                                                         text = element_text(size=30),
                                                         strip.text.x = element_text(size = 25),
                                                         strip.text.y = element_text(size = 25),
                                                         legend.text=element_text(size = 25)) +
  facet_grid(SelfReg~SigmaR, scales = "free",  labeller = label_bquote(cols = sigma[R] == .(SigmaR), rows = .(SelfReg))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[B])),
       y = expression("Coexisting Fraction"~(phi)),
       color = expression("Mean\nInteraction\nStrength"(mu[B])),
       shape = expression("Mean\nInteraction\nStrength"(mu[B])))
plSICubic

jpeg("../CavityHOIs-Paper/figs/SIFigCubic.jpeg", width = 4000, height = 2100, res = 300)
plSICubic
dev.off()

