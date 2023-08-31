library(ggplot2)
library(fitdistrplus)
library(sjPlot)
library(ggpubr)
library(plyr)
library(dplyr)
library(tidyr)

# Change margin settings
par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/paml/")

introgression_tract_dN_90 = "mean_dN_introgression_tract_90.csv"
species_tree_tract_dN_90 = "mean_dN_species_tree_90.csv"

introgression_tract_dS_90 = "mean_dS_introgression_tract_90.csv"
species_tree_tract_dS_90 = "mean_dS_species_tree_90.csv"

introgression_tract_dNdS_90 = "mean_dNdS_introgression_tract_90.csv"
species_tree_tract_dNdS_90 = "mean_dNdS_species_tree_90.csv"

# Read the CSV file as a one-dimensional vector
dist_dN_introgression_tract_90 <- scan(introgression_tract_dN_90, what = numeric(), sep = ",")
dist_dN_species_tree_90 <- scan(species_tree_tract_dN_90, what = numeric(), sep = ",")
dist_dS_introgression_tract_90 <- scan(introgression_tract_dS_90, what = numeric(), sep = ",")
dist_dS_species_tree_90 <- scan(species_tree_tract_dS_90, what = numeric(), sep = ",")
dist_dNdS_introgression_tract_90 <- scan(introgression_tract_dNdS_90, what = numeric(), sep = ",")
dist_dNdS_species_tree_90 <- scan(species_tree_tract_dNdS_90, what = numeric(), sep = ",")

df <- data.frame(
  dist_type = factor(c(rep("introgression_tract", 1000), rep("species_tree_tract", 1000))),
  mean_dN = c(dist_dN_introgression_tract_90, dist_dN_species_tree_90),
  mean_dS = c(dist_dS_introgression_tract_90, dist_dS_species_tree_90),
  mean_dNdS = c(dist_dNdS_introgression_tract_90, dist_dNdS_species_tree_90))

# Convert data to long format
df_long <- df %>% pivot_longer(cols = c(mean_dN, mean_dS, mean_dNdS),
                               names_to = "metric",
                               values_to = "value")

library(ggplot2)

# Your data
df <- data.frame(
  dist_type = factor(c(rep("introgression_tract", 1000), rep("species_tree_tract", 1000))),
  mean_dN = c(dist_dN_introgression_tract_90, dist_dN_species_tree_90),
  mean_dS = c(dist_dS_introgression_tract_90, dist_dS_species_tree_90),
  mean_dNdS = c(dist_dNdS_introgression_tract_90, dist_dNdS_species_tree_90))

# Reorder the levels of the metric variable
df_long$metric <- factor(df_long$metric, levels = c("mean_dN", "mean_dS", "mean_dNdS"))

# Create the density plots
fig1 <- ggplot(df_long, aes(x = value, fill = dist_type)) +
  geom_density(aes(y = after_stat(density)), alpha = 0.8) +
  facet_wrap(
    ~ metric, 
    scales = "free", 
    strip.position = "bottom",
    labeller = labeller(
      metric = c(
        mean_dN = "Mean dN",
        mean_dS = "Mean dS",
        mean_dNdS = "Mean dNdS"
      )
    )
  ) +
  scale_fill_manual(values = c("introgression_tract" = "#A6CEE3", "species_tree_tract" = "#FDBF6F"),
                    labels = c("Introgression Tract", "Genome-Wide Background")) + 
  labs(x = "", y = "Probability Density", fill = "") +
  theme_blank() +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "right",
    panel.grid = element_blank()
  ) + 
  guides(fill = guide_legend(title = NULL))

# Print the plot
fig1


ggsave(filename = "dist_paml.eps", plot = fig1)
ggsave(filename = "dist_paml.png", plot = fig1)

mean(dist_dN_introgression_tract_90)
mean(dist_dN_species_tree_90)
mean(dist_dS_introgression_tract_90)
mean(dist_dS_species_tree_90)
mean(dist_dNdS_introgression_tract_90)
mean(dist_dNdS_species_tree_90)

####
T-Tests
t.test(mean_dN~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
t.test(mean_dS~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
t.test(mean_dNdS~dist_type, data=df, paired=FALSE, var.eq=F, alternative="two.sided")
