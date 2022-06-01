
source("Functions.R")

# reading in the data
out_abds <- read.csv("simdata/sim0110_Mean.csv", header = TRUE)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)
abd_stats <- GetStatistics(proc_abds)

melt_stats <- abd_stats %>%
  ungroup() %>%
  mutate(SelfReg = ifelse(MuD == 1, "Constant", "Mean-Adjusted")) %>%
  dplyr::select(S, SelfReg, Mu, Sigma, Phi) %>%
  melt(id.vars = c("S", "SelfReg", "Mu", "Sigma")) %>%
  mutate(value = ifelse(Mu == min(Mu), -value, value)) %>%
  dplyr::group_by(S, SelfReg) %>%
  dplyr::summarise(Difference = sum(value))

plSIMean <- ggplot(melt_stats, aes(x = S, y = Difference, color = SelfReg, shape = SelfReg)) +
  geom_point(size = 7, alpha = 0.9) + theme_bw() +
  geom_hline(yintercept = 0, color = "red") +
  theme(legend.position = "right",
        panel.spacing = unit(1.5, "lines"),
        text = element_text(size=20),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        legend.text=element_text(size = 15)) +
  labs(x = expression("Number of Species"~(S)),
       y = expression(atop(" Difference in Fractions of","Coexisting Species"~(phi))),
       color = expression("Self-Regulation"~ (mu[D])),
       shape = expression("Self-Regulation"~ (mu[D])))
plSIMean

jpeg("../CavityHOIs-Paper/figs/SIFigMean.jpeg", width = 2500, height = 1500, res = 300)
plSIMean
dev.off()

