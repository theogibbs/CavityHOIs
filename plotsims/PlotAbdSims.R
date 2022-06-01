source("Functions.R")

S <- 50

mu_R <- 1.5
sigma_R <- 0

r <- rnorm(S, mean = mu_R, sd = sigma_R)

A <- BuildA(S = S, mu_A = 0, sigma_A = 1, rho_A = 0)
diag(A) <- 0
A[A != 0] <- A[A != 0] / sd(A[A != 0]) / sqrt(S)
A[A != 0] <- A[A != 0] - mean(A[A != 0])

B <- array(rnorm(S*S*S, mean = 0, sd = 1), c(S, S, S))
B <- NoSelfHOIs(B)

B[B != 0] <- B[B != 0] / sd(B[B != 0]) / S
B[B != 0] <- B[B != 0] - mean(B[B != 0])

mu_len <- 3
sigma_len <- 50

mu_As <- seq(-1, -2, length.out = mu_len)
sigma_As <- seq(0.01, 0.85, length.out = sigma_len) # Code breaks if sigma = 0

mu_Bs <- seq(-1, -2, length.out = mu_len)
sigma_Bs <- seq(0.01, 0.85, length.out = sigma_len)  # Code breaks if sigma = 0

settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

out_data <- data.frame(Mu = c(), Sigma = c(), Interaction = c(), SpeciesID = c(),
                       Abundances = c(), Equilibrium = c(), Uninvadeable = c())

for(mu in 1:mu_len) {
  print(mu)
  for(sigma in 1:sigma_len) {
    print(sigma)
    
    curB <- B
    curB[curB != 0] <- curB[curB != 0] * sigma_Bs[sigma]
    curB[curB != 0] <- curB[curB != 0] + mu_Bs[mu] / S^2
    
    pars <- list(S = S,
                 r = r,
                 d = rep(1, times = S),
                 A = matrix(0, nrow = S, ncol = S),
                 B = curB,
                 h = 0)
    cur_hoi_abds <- GetAbds(pars, settings)
    
    curN <- cur_hoi_abds$Abundances
    curN[curN > 0] <- 1
    curN[curN < 0.01] <- 0
    
    cur_ints <- curB
    
    for(i in 1:S) {
      for(j in 1:S) {
        for(k in 1:S) {
          cur_ints[i, j, k] <- curN[i] * curN[j] * curN[k] * cur_ints[i, j, k]
        }
      }
    }
    
    cur_ints[cur_ints == 0] <- NA
    
    cur_mean <- S^2 * mean(cur_ints, na.rm = TRUE)
    cur_sd <- S * sd(cur_ints, na.rm = TRUE)
    
    hoi_data <- cbind(data.frame(Mu = mu_Bs[mu], Sigma = sigma_Bs[sigma],
                                 MeanInt = cur_mean, SdInt = cur_sd,
                                 Interaction = "Higher Order"), cur_hoi_abds)
    
    curA <- A
    curA[curA != 0] <- curA[curA != 0] * sigma_As[sigma]
    curA[curA != 0] <- curA[curA != 0] + mu_As[mu] / S
    
    pars <- list(S = S,
                 r = r,
                 d = rep(1, times = S),
                 A = curA,
                 B = array(0, c(S, S, S)),
                 h = 0)
    cur_pw_abds <- GetAbds(pars, settings)
    
    curN <- cur_pw_abds$Abundances
    curN[curN > 0] <- 1
    curN[curN < 0.01] <- 0
    curN <- as.logical(outer(curN, curN))
    cur_ints <- curA[curN]
    cur_ints[cur_ints == 0] <- NA
    cur_mean <- S * mean(cur_ints, na.rm = TRUE)
    cur_sd <- sqrt(S) * sd(cur_ints, na.rm = TRUE)
    
    pw_data <- cbind(data.frame(Mu = mu_As[mu], Sigma = sigma_As[sigma],
                                MeanInt = cur_mean, SdInt = cur_sd,
                                Interaction = "Pairwise"), cur_pw_abds)
    
    out_data <- rbind(out_data, pw_data, hoi_data)
    
  }
}

plot_data <- out_data %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order"))) %>%
  mutate(SpeciesID = as.factor(SpeciesID)) %>%
  mutate(Mu = - Mu)

plAbds <- ggplot(plot_data, aes(x = Sigma, y = Abundances, color = SpeciesID)) +
  geom_line(size = 0.75, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = "none",
        panel.spacing = unit(1, "lines"),
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  geom_abline(slope = 0, intercept = mu_R) +
  geom_abline(slope = 0, intercept = 0, color = "red") +
  facet_grid(Interaction ~ Mu, labeller = label_bquote(cols = mu[A]  ~ "or" ~ mu[B] == .(Mu), rows = .(Interaction))) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A] ~ "or" ~ sigma[B])),
       y = "Abundance")
plAbds

jpeg("../CavityHOIs-Paper/figs/SIFigAbundances.jpeg", width = 2000, height = 1250, res = 300)
plAbds
dev.off()

plot_data <- out_data %>%
  select(Mu, Sigma, Interaction, MeanInt, SdInt) %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order"))) %>%
  mutate(Mu = as.factor(- Mu)) %>%
  mutate(MeanInt = - MeanInt) %>%
  melt(id.vars = c("Mu", "Sigma", "Interaction")) %>%
  mutate(variable = ifelse(variable == "MeanInt", "Realized Mean", "Realized Standard Deviation"))

plAbdInts <- ggplot(plot_data, aes(x = Sigma, y = value, color = Mu)) +
  geom_line(size = 0.75, alpha = 0.5) +
  theme_bw() +
  theme(panel.spacing = unit(1, "lines"),
        text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  facet_wrap(Interaction ~ variable, scales = "free") +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A] ~ "or" ~ sigma[B])),
       y = "",
       color = expression(atop("Mean\nInteraction\nStrength",(mu[A] ~"or" ~mu[B]))))
plAbdInts

jpeg("../CavityHOIs-Paper/figs/SIFigAbdInts.jpeg", width = 2750, height = 1500, res = 300)
plAbdInts
dev.off()
