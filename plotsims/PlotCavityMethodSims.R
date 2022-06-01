
source("Functions.R")

# reading in the data

out_abds <- read.csv("simdata/sim1204_CavityMethod.csv", header = TRUE)

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
  mutate(Mu = as.factor(-Mu)) %>%
  filter(Interaction == "Higher Order")

plCavityMethod <- ggplot(melt_stats, aes(x = Sigma, y = value, color = Mu, shape = Mu)) +
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
       shape = expression(atop("Mean\nInteraction\nStrength", (mu[B]))))
plCavityMethod

jpeg("../CavityHOIs-Paper/figs/FigCavityMethod.jpeg", width = 4500, height = 1500, res = 300)
plCavityMethod
dev.off()

# More diverse ecosystems

fileregex <- "sim522_PairCorr."
out_abds_corr <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))


fileregex <- "sim51._AllInts300."
out_abds_all <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE)) %>%
  mutate(ParsID = 0)

out_abds <- rbind(out_abds_all, out_abds_corr)

CheckAbds(out_abds)
out_abds <- out_abds %>%
  filter(Equilibrium == TRUE, Uninvadeable == TRUE) %>%
  mutate(h = 0)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)
pred_stats <- GetPredictions(abd_stats, 10)

melt_stats <- pred_stats %>%
  ungroup() %>%
  dplyr::select(S, Mu, Sigma, RhoA, Phi, MeanAbd, VarAbd, Interaction) %>%
  melt(id.vars = c("S", "Mu", "Sigma", "RhoA", "Interaction")) %>%
  mutate(Error = c(pred_stats$ErrorPhi, pred_stats$ErrorMean, pred_stats$ErrorVar)) %>%
  mutate(Prediction = c(pred_stats$PredFraction, pred_stats$PredMean, pred_stats$PredSec - pred_stats$PredMean^2)) %>%
  mutate(Mu = as.factor(-Mu)) %>%
  mutate(RhoA = as.factor(RhoA)) %>%
  mutate(variable = c(rep("Coexisting Fraction", times = nrow(pred_stats)),
                      rep("SAD Mean", times = nrow(pred_stats)),
                      rep("SAD Variance", times = nrow(pred_stats))))

plDivCM <- ggplot(melt_stats, aes(x = Sigma, y = value, color = RhoA, shape = RhoA)) +
  geom_errorbar(aes(ymin = value - Error, ymax = value + Error), width = 0, size = 2) +
  geom_point(size = 7) + theme_bw() + theme(legend.position = "top",
                                            panel.spacing = unit(1.5, "lines"),
                                            text = element_text(size=30),
                                            strip.text.x = element_text(size = 25),
                                            strip.text.y = element_text(size = 20),
                                            legend.text=element_text(size = 25)) +
  geom_line(aes(x = Sigma, y = Prediction, color = RhoA), size = 3, alpha = 0.6) +
  facet_wrap(Interaction~variable, scales = "free") +
  labs(x = expression("Variation in Interaction Strengths"~(sqrt(sigma[A]^2+sigma[B]^2))),
       y = " ",
       color = expression("Pairwise Correlation"~(rho[A])),
       shape = expression("Pairwise Correlation"~(rho[A])))
plDivCM

jpeg("../CavityHOIs-Paper/figs/SIFigDivCM.jpeg", width = 4500, height = 3750, res = 300)
plDivCM
dev.off()



