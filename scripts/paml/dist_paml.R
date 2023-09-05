library(ggplot2)
library(fitdistrplus)
library(sjPlot)
library(ggpubr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggh4x)

# Change margin settings
#par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/paml/")

introgression_tract_dN_90 = "mean_dN_introgression_tract_90.csv"
species_tree_tract_dN_90 = "mean_dN_species_tree_90.csv"

introgression_tract_dS_90 = "mean_dS_introgression_tract_90.csv"
species_tree_tract_dS_90 = "mean_dS_species_tree_90.csv"

introgression_tract_dNdS_90 = "mean_dNdS_introgression_tract_90.csv"
species_tree_tract_dNdS_90 = "mean_dNdS_species_tree_90.csv"

introgression_tract_dN_80 = "mean_dN_introgression_tract_80.csv"
species_tree_tract_dN_80 = "mean_dN_species_tree_80.csv"

introgression_tract_dS_80 = "mean_dS_introgression_tract_80.csv"
species_tree_tract_dS_80 = "mean_dS_species_tree_80.csv"

introgression_tract_dNdS_80 = "mean_dNdS_introgression_tract_80.csv"
species_tree_tract_dNdS_80 = "mean_dNdS_species_tree_80.csv"

# Read the CSV file as a one-dimensional vector
dist_dN_introgression_tract_90 <- scan(introgression_tract_dN_90, what = numeric(), sep = ",")
dist_dN_species_tree_90 <- scan(species_tree_tract_dN_90, what = numeric(), sep = ",")
dist_dS_introgression_tract_90 <- scan(introgression_tract_dS_90, what = numeric(), sep = ",")
dist_dS_species_tree_90 <- scan(species_tree_tract_dS_90, what = numeric(), sep = ",")
dist_dNdS_introgression_tract_90 <- scan(introgression_tract_dNdS_90, what = numeric(), sep = ",")
dist_dNdS_species_tree_90 <- scan(species_tree_tract_dNdS_90, what = numeric(), sep = ",")

dist_dN_introgression_tract_80 <- scan(introgression_tract_dN_80, what = numeric(), sep = ",")
dist_dN_species_tree_80 <- scan(species_tree_tract_dN_80, what = numeric(), sep = ",")
dist_dS_introgression_tract_80 <- scan(introgression_tract_dS_80, what = numeric(), sep = ",")
dist_dS_species_tree_80 <- scan(species_tree_tract_dS_80, what = numeric(), sep = ",")
dist_dNdS_introgression_tract_80 <- scan(introgression_tract_dNdS_80, what = numeric(), sep = ",")
dist_dNdS_species_tree_80 <- scan(species_tree_tract_dNdS_80, what = numeric(), sep = ",")

df <- data.frame(
  dist_type = factor(c(rep("introgression_tract", 1000), rep("species_tree_tract", 1000), rep("introgression_tract", 1000), rep("species_tree_tract", 1000))),
  mean_dN = c(dist_dN_introgression_tract_90, dist_dN_species_tree_90, dist_dN_introgression_tract_80, dist_dN_species_tree_80),
  mean_dS = c(dist_dS_introgression_tract_90, dist_dS_species_tree_90, dist_dS_introgression_tract_80, dist_dS_species_tree_80),
  mean_dNdS = c(dist_dNdS_introgression_tract_90, dist_dNdS_species_tree_90, dist_dNdS_introgression_tract_80, dist_dNdS_species_tree_80),
  posterior_probability = factor(c(rep("probability_90", 2000), rep("probability_80", 2000)))
  )
  
# Convert data to long format
df_long <- df %>% pivot_longer(cols = c(mean_dN, mean_dS, mean_dNdS),
                               names_to = "metric",
                               values_to = "value")

# Reorder the levels of the metric variable
df_long$metric <- factor(df_long$metric, levels = c("mean_dN", "mean_dS", "mean_dNdS"))

# Reorder the levels of the metric variable
df_long$posterior_probability <- factor(df_long$posterior_probability, levels = c("probability_90", "probability_80"))

###########

# Create the ggplot using facet_grid for independent scales
fig1 <- ggplot(df_long, aes(x = value, fill = dist_type)) +
  geom_density(aes(y=after_stat(density)), alpha=0.8) +
  facet_grid2(
    posterior_probability ~ metric,
    scales = "free", independent = "all",
    labeller = labeller(posterior_probability = c(
      "probability_90" = "Posterior\nProbability\n> 90%",
      "probability_80" = "Posterior\nProbability\n> 80%"),
      metric = c(
        mean_dN = "Mean dN",
        mean_dS = "Mean dS",
        mean_dNdS = "Mean dNdS"
      )
    ), switch = "x") +
  scale_fill_manual(values = c("introgression_tract" = "#A6CEE3", "species_tree_tract" = "#FDBF6F"),
                    labels = c("Introgression Tract", "Genome-Wide Background")) +
  labs(x = "", y = "Probability Density", fill = "") +
  theme_blank() +
  theme(strip.placement = "outside",
        axis.title = element_text(size = 20),
        legend.position = "top",
        panel.grid = element_blank(),
        strip.text.y = element_text(size = 12, angle = 0),
        strip.text.x = element_text(size = 16), # Font size for facet labels
        legend.text = element_text(size = 16)
  ) +
  guides(fill = guide_legend(title = NULL))

fig1

ggsave(filename = "dist_paml.eps", plot = fig1)
ggsave(filename = "dist_paml.png", plot = fig1)

mean(dist_dN_introgression_tract_90)
mean(dist_dN_species_tree_90)
mean(dist_dS_introgression_tract_90)
mean(dist_dS_species_tree_90)
mean(dist_dNdS_introgression_tract_90)
mean(dist_dNdS_species_tree_90)

mean(dist_dN_introgression_tract_80)
mean(dist_dN_species_tree_80)
mean(dist_dS_introgression_tract_80)
mean(dist_dS_species_tree_80)
mean(dist_dNdS_introgression_tract_80)
mean(dist_dNdS_species_tree_80)

####
T-Tests
t.test(mean_dN~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
t.test(mean_dS~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
t.test(mean_dNdS~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
