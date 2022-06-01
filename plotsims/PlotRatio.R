sigma_R <- 0
sigma_A <- 0.1
sigma_B <- 0.1

mu_0s <- seq(0.001, 1, length.out = 100)

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)

plot_data <- data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R)

sigma_R <- 0.5

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))

sigma_R <- 0.25

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))

sigma_R <- 0.05

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))


sigma_R <- 0.01

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))


sigma_R <- 0.1

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))


sigma_R <- 1

ratios <- mu_0s / sqrt((1 - sigma_A^2 - sqrt((1-sigma_A^2)^2 - 4 * sigma_B^2 * (sigma_R^2 + mu_0s^2))) / (2 * sigma_B^2) - mu_0s^2)
plot_data <- rbind(plot_data, data.frame(Mu0 = mu_0s, Ratio = ratios, SigmaR = sigma_R))


plot_data <- plot_data %>%
  mutate(SigmaR = as.factor(SigmaR))

ggplot(plot_data, aes(x = Mu0, y = Ratio, color = SigmaR)) + 
  geom_line() + theme_bw()



plRatio <- ggplot(plot_data, aes(x = Mu0, y = Ratio, color = SigmaR)) +
  geom_line(size = 3, alpha = 0.8) + theme_bw() + theme(legend.position = "right",
                                           text = element_text(size=30),
                                           strip.text.x = element_text(size = 25),
                                           strip.text.y = element_text(size = 20),
                                           legend.text=element_text(size = 25)) +
  labs(x = expression("Mean Abundance"~(mu[0])),
       y = expression("Ratio"~(mu[0] / sigma[0])),
       color = expression(sigma[R]))
plRatio


jpeg("../CavityHOIs-Paper/figs/SIFigRatio.jpeg", width = 2500, height = 2250, res = 300)
plRatio
dev.off()

