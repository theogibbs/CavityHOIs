

out_sigmas <- read.csv("simdata/sim1204_MayBound.csv", header = TRUE)

out_sigmas %>% filter(MaxSteps == TRUE)

# make sure this is consistent with the parameters used to run the simulation
max_sigma <- 1
num_sigmas <- 2e3
step_size <- max_sigma / num_sigmas
print(step_size)
test_sigmas <- seq(0, max_sigma, length.out = num_sigmas)

out_sigmas %>% filter(Sigma < test_sigmas[3])

Ss <- unique(out_sigmas$S)
num_repl <- 1e3

maxs <- c()
for(S in Ss) {
  cur_maxs <- c()
  for(i in 1:num_repl) {
    cur_maxs <- c(cur_maxs, max(rnorm(S)))
  }
  maxs <- c(maxs, mean(cur_maxs))
}

plot(Ss, maxs)
points(Ss,  sqrt(2 * log(Ss)))

mean_sigmas <- out_sigmas %>%
  dplyr::group_by(S, MuR, SigmaR, MuD, SigmaD, MuA, SigmaA, RhoA, MuB, SigmaB, RhoB, h, Interaction) %>%
  dplyr::summarise(MinSigma = min(Sigma), MaxSigma = max(Sigma), NumRepl = length(Sigma), Sigma = mean(Sigma)) %>%
  dplyr::mutate(Const = maxs[which(S == Ss)]) %>%
  dplyr::mutate(Prediction = ifelse(Interaction == "Pairwise", sqrt(1 / S / (1 + Const^2)), 1 / S * Const / (1 + Const^2) / MuR)) %>%
  dplyr::mutate(Interaction = ifelse(Interaction == "Higher Order", "Higher-order", Interaction)) %>%
  dplyr::mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher-order")))

print(unique(mean_sigmas$NumRepl))

plMayBound <- ggplot(mean_sigmas, aes(x = S, y = Sigma, color = Interaction, shape = Interaction)) +
  geom_errorbar(aes(ymin = MinSigma, ymax = MaxSigma), alpha = 0.6, width = 0, size = 2) +
  geom_line(aes(x = S, y = Prediction, color = Interaction), size = 2, alpha = 0.6) +
  geom_point(size = 7, alpha = 0.8) + theme_classic() +
  theme(legend.position = c(0.8, 0.85),
        text = element_text(size=30),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size = 25)) +
  labs(x = expression("Number of Species"~(S)),
       y = expression("Critical Variability"~(tilde(sigma)[A]~"or"~tilde(sigma)[B])),
       color = "Interaction",
       shape = "Interaction") + scale_y_log10()
plMayBound

jpeg("../CavityHOIs-Paper/figs/FigMayBound.jpeg", width = 2500, height = 2250, res = 300)
plMayBound
dev.off()

