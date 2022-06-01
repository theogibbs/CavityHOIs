
source("Functions.R")

# reading in the data
out_abds <- read.csv("simdata/sim1215_MultAtt.csv", header = TRUE)
new_abds <- read.csv("simdata/sim_SatMultAtt.csv", header = TRUE)
out_abds <- rbind(out_abds, new_abds)

att_data <- LabelAbds(out_abds) %>%
  dplyr::group_by(Sigma, h, Interaction, Mu, ParsID, SpeciesID) %>%
  dplyr::mutate(StDevAbd = sd(Abundances)) %>%
  dplyr::group_by(Sigma, h, Interaction, Mu, ParsID) %>%
  dplyr::summarise(SumStDev = mean(StDevAbd), UFP = (Equilibrium * Uninvadeable)) %>%
  dplyr::mutate(SumStDev = ifelse(UFP, SumStDev, NA)) %>%
  dplyr::mutate(SumStDev = ifelse(SumStDev < 1e10, SumStDev, NA)) %>%
  dplyr::group_by(Sigma, h, Interaction, Mu) %>%
  dplyr::summarise(SumStDev = mean(SumStDev, na.rm = TRUE), UFP = mean(UFP)) %>%
  melt(id.vars = c("Mu", "h", "Sigma", "Interaction")) %>%
  mutate(variable = ifelse(variable == "SumStDev", "(A) Average Standard Deviation",
                           "(B) Probability of Non-invadeable Equilibrium"))

plMultAtt <- ggplot(att_data %>% filter(h == 0),
                    aes(x = Sigma, y = value, color = Interaction, shape = Interaction)) +
  facet_wrap(~variable, scales = "free") +  theme_bw() +
  theme(legend.position = "right",
        panel.spacing = unit(1, "lines"),
        text = element_text(size=25),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        legend.text=element_text(size = 25)) +
  geom_point(size = 6) +
  geom_line(size = 2, alpha = 0.5) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A]~"or"~sigma[B])), y = " ")
plMultAtt

jpeg("../CavityHOIs-Paper/figs/SIFigMultAtt.jpeg", width = 4000, height = 1500, res = 300)
plMultAtt
dev.off()

plSatMultAtt <- ggplot(att_data %>% filter(h != 0), aes(x = Sigma, y = value)) +
  facet_wrap(~variable, scales = "free") +  theme_bw() +
  theme(legend.position = "right",
        panel.spacing = unit(1, "lines"),
        text = element_text(size=25),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 20),
        legend.text=element_text(size = 25)) +
  geom_point(size = 6) +
  geom_line(size = 2, alpha = 0.5) +
  labs(x = expression("Variation in Interaction Strengths"~(sigma[A]~"or"~sigma[B])), y = " ")
plSatMultAtt

jpeg("../CavityHOIs-Paper/figs/SIFigSatMultAtt.jpeg", width = 4000, height = 1500, res = 300)
plSatMultAtt
dev.off()

