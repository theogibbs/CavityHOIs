source("Functions.R")

### making the dynamics plot

input_S <- 20
input_mu_r <- 1.5
input_sigma_r <- 0.5
input_mu_d <- 1
input_sigma_d <- 0
input_mu_A <- -1
input_sigma_A <- 0.5
input_rho_A <- 0
input_mu_B <- -1
input_sigma_B <- 0.25
input_rho_B <- 0
input_h <- 0

input_params <- crossing(S = input_S, MuR = input_mu_r,
                         SigmaR = input_sigma_r, MuD = input_mu_d,
                         SigmaD = input_sigma_d, MuA = input_mu_A,
                         SigmaA = input_sigma_A, RhoA = input_rho_A,
                         MuB = input_mu_B, SigmaB = input_sigma_B,
                         RhoB = input_rho_B, h = input_h,
                         ParsID = 1, CommunityID = 1)

pars <- BuildPars(input_params)

ini_state <- runif(pars$S, min = 0, max = 0.01)
end_time <- 20
time_step <- 0.01
out_dyn <- IntegrateDynamics(ini_state, pars, end_time, time_step, Dynamics)
plSeries <- PlotSeries(out_dyn, title = "(B) Community dynamics")
plSeries


jpeg("../CavityHOIs-Paper/figs/FigDynamics.jpeg", width = 2500, height = 2000, res = 300)
plSeries
dev.off()

# making histogram plot
fileregex <- "sim508_OnlyPairwise300."
out_abds_pw <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

fileregex <- "sim510_OnlyHOIs300."
out_abds_hois <- paste0("simdata/", list.files(path = "./simdata/", pattern = fileregex)) %>%
  map_df(~read.csv(., header = TRUE))

out_abds <- rbind(out_abds_pw, out_abds_hois)

CheckAbds(out_abds)
out_abds <- out_abds %>% filter(Equilibrium == TRUE, Uninvadeable == TRUE)
proc_abds <- LabelAbds(out_abds)

hist_sigmas <- c(0.5, 1)

plot_abds <- proc_abds %>%
  filter(Sigma == hist_sigmas[1] | Sigma == hist_sigmas[2]) %>%
  filter(Abundances > 0) %>%
  filter(MuB == 0 & SigmaB == 0 | MuA == 0 & SigmaA == 0)

pred_data <- plot_abds %>%
  dplyr::select(-one_of("Abundances", "SpeciesID", "CommunityID", "Index")) %>%
  unique()

pred_abds <- GetPredictions(pred_data)
pred_abds <- GetAdjPreds(pred_abds)

plot_abds <- merge(plot_abds, pred_abds)
plot_abds <- plot_abds %>%
  mutate(PredDensity = dtruncnorm(Abundances, a = 0, b = Inf, mean = AdjPredMean, sd = AdjPredSD))

plHist <- ggplot(plot_abds, aes(x = Abundances, y = ..density..)) +
  geom_histogram(fill = "white", color = "black", binwidth = 0.05) +
  ggtitle("(C) Species abundance distributions") +
  facet_grid(Sigma ~ Interaction, labeller = label_bquote(cols = .(Interaction), rows = sigma[A]~ "or"~ sigma[B] == .(Sigma))) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 25),
                     strip.text.y = element_text(size = 25),
                     text = element_text(size=30)) + ylab("Density") +
  geom_line(aes(x = Abundances, y = PredDensity), color = "blue", size = 1) #+
  #geom_vline(data = filter(plot_abds, Sigma == 1, Interaction == "Pairwise"),
  #           aes(xintercept = filter(pred_abds, Sigma == 1, Interaction == "Pairwise")$AdjPredMean), color = "darkgreen") +
  #geom_vline(data = filter(plot_abds, Sigma == 1, Interaction == "Pairwise"),
  #           aes(xintercept = filter(pred_abds, Sigma == 1, Interaction == "Pairwise")$AdjPredMean +
  #                 filter(pred_abds, Sigma == 1, Interaction == "Pairwise")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Pairwise"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Pairwise")$AdjPredMean), color = "darkgreen") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Pairwise"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Pairwise")$AdjPredMean +
  #                 filter(pred_abds, Sigma == 0.5, Interaction == "Pairwise")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Pairwise"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Pairwise")$AdjPredMean -
  #                 filter(pred_abds, Sigma == 0.5, Interaction == "Pairwise")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Higher Order")$AdjPredMean), color = "darkgreen") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Higher Order")$AdjPredMean +
  #                 filter(pred_abds, Sigma == 0.5, Interaction == "Higher Order")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 0.5, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 0.5, Interaction == "Higher Order")$AdjPredMean -
  #                 filter(pred_abds, Sigma == 0.5, Interaction == "Higher Order")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 1, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 1, Interaction == "Higher Order")$AdjPredMean), color = "darkgreen") +
  #geom_vline(data = filter(plot_abds, Sigma == 1, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 1, Interaction == "Higher Order")$AdjPredMean +
  #                 filter(pred_abds, Sigma == 1, Interaction == "Higher Order")$AdjPredSD), color = "darkgreen", linetype = "dashed") +
  #geom_vline(data = filter(plot_abds, Sigma == 1, Interaction == "Higher Order"),
  #           aes(xintercept = filter(pred_abds, Sigma == 1, Interaction == "Higher Order")$AdjPredMean -
  #                 filter(pred_abds, Sigma == 1, Interaction == "Higher Order")$AdjPredSD), color = "darkgreen", linetype = "dashed")
plHist

jpeg("../CavityHOIs-Paper/figs/FigSAD.jpeg", width = 3000, height = 2000, res = 300)
plHist
dev.off()

