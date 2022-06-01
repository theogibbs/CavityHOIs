source("Functions.R")


settings <- list(abd_cutoff = 1e-14, gr_cutoff = 0.01, endtime = 1e7, inimin = 0, inimax = 1)

# some scratch space and code to generate the desired combination of parameters
input_S <- 200
input_mu_r <- 1.5
input_sigma_r <- 0
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- 0.25
input_rho_A <- 0
input_mu_B <- 0
input_sigma_B <- 0
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h, ParsID = 1)

pw_pars <- BuildPars(input_params)
curA <- pw_pars$A
pw_abds <- GetAbds(pw_pars, settings) %>%
  mutate(TotInts = rowSums(curA)) %>%
  mutate(Interaction = "Pairwise")

# higher-order interactions
input_mu_A <- 0
input_sigma_A <- 0
input_mu_B <- -1
input_sigma_B <- 0.35

input_params <- crossing(S = input_S, MuR = input_mu_r, SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A, SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B, RhoB = input_rho_B, h = input_h, ParsID = 1)

hoi_pars <- BuildPars(input_params)
curB <- hoi_pars$B
hoi_abds <- GetAbds(hoi_pars, settings) %>%
  mutate(TotInts = rowSums(curB)) %>%
  mutate(Interaction = "Higher Order")

abd_data <- rbind(pw_abds, hoi_abds) %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order")))

plTotInts <- ggplot(abd_data, aes(x = Abundances, y = TotInts)) +
  theme_bw() + stat_smooth(formula = y ~ x, method = "lm") +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  geom_point(size = 2, alpha = 0.5) + facet_wrap(~Interaction, scales = "free") +
  labs(x = "Abundance", y = "Focal Species\nTotal Interaction Strength")
plTotInts

jpeg("../CavityHOIs-Paper/figs/SIFigTotalInts.jpeg", width = 2250, height = 1250, res = 300)
plTotInts
dev.off()

pw_abds <- pw_abds$Abundances
pw_ints <- melt(-curA) %>%
  mutate(FocAbd = pw_abds[Var1]) %>%
  mutate(IntAbd = pw_abds[Var2]) %>%
  select(FocAbd, IntAbd, value) %>%
  mutate(Interaction = "Pairwise")

hoi_abds <- hoi_abds$Abundances
hoi_ints <- melt(-curB) %>%
  mutate(FocAbd = hoi_abds[Var1]) %>%
  mutate(IntAbd = hoi_abds[Var2] * hoi_abds[Var3]) %>%
  select(FocAbd, IntAbd, value) %>%
  mutate(Interaction = "Higher Order")

int_data <- rbind(pw_ints, hoi_ints) %>%
  mutate(Interaction = factor(Interaction, levels = c("Pairwise", "Higher Order")))

low_ints <- int_data %>%
  group_by(Interaction) %>%
  filter(FocAbd == min(FocAbd))

plLow <- ggplot(low_ints, aes(x = IntAbd, y = value)) +
  theme_bw() + stat_smooth(formula = y ~ x, method = "lm") +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  geom_point(size = 2, alpha = 0.1) + facet_wrap(~Interaction, scales = "free") +
  labs(x = "Interacting Abundance Strength", y = "Interaction Strength") +
  ggtitle("(A) Low abundance species")
plLow

high_ints <- int_data %>%
  group_by(Interaction) %>%
  filter(FocAbd == max(FocAbd))

plHigh <- ggplot(high_ints, aes(x = IntAbd, y = value)) +
  theme_bw() + stat_smooth(formula = y ~ x, method = "lm") +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  geom_point(size = 2, alpha = 0.1) + facet_wrap(~Interaction, scales = "free") +
  labs(x = "Interacting Abundance Strength", y = "Interaction Strength") +
  ggtitle("(B) High abundance species")
plHigh

slope_data <- data.frame(Interaction = c(), FocAbd = c(), Slope = c())
for(i in 1:length(pw_abds)) {
  cur_pw <- pw_ints %>%
    filter(FocAbd == pw_abds[i])
  pw_slope <- summary(lm(value ~ IntAbd, data = cur_pw))$coefficients[2,1]
  pw_slope <- data.frame(Interaction = "Pairwise",
                         FocAbd = pw_abds[i],
                         Slope = pw_slope)
  
  cur_hoi <- hoi_ints %>%
    filter(FocAbd == hoi_abds[i])
  hoi_slope <- summary(lm(value ~ IntAbd, data = cur_hoi))$coefficients[2,1]
  hoi_slope <- data.frame(Interaction = "Higher Order",
                          FocAbd = hoi_abds[i],
                          Slope = hoi_slope)
  
  slope_data <- rbind(slope_data, pw_slope, hoi_slope)
}

plSlopes <- ggplot(slope_data, aes(x = FocAbd, y = Slope)) +
  theme_bw() + stat_smooth(formula = y ~ x, method = "lm") +
  theme(text = element_text(size=20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text=element_text(size=15)) +
  geom_point(size = 2, alpha = 0.5) + facet_wrap(~Interaction, scales = "free") +
  labs(x = "Focal Abundance", y = "Slope") +
  ggtitle("(C) Regressions for all abundances")
plSlopes

jpeg("../CavityHOIs-Paper/figs/SIFigInteractions.jpeg", width = 3500, height = 2500, res = 300)
grid.arrange(plLow, plHigh, plSlopes, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()
