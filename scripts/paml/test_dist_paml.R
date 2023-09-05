library(ggplot2)
library(tidyr)
library(gridExtra)

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

library(ggplot2)
library(tidyr)

# Reshape the data to long format
df_long <- df %>%
  gather(variable, value, mean_dN, mean_dS, mean_dNdS)

# Create a ggplot2 density plot
density_plot <- ggplot(df_long, aes(x = value, fill = dist_type)) +
  geom_density(alpha = 0.7) +
  labs(title = "Density Plots",
       x = "Value",
       y = "Density") +
  facet_grid(posterior_probability ~ variable, scales = "free_x") + 
  scale_fill_manual(values = c("introgression_tract" = "#A6CEE3", "species_tree_tract" = "#FDBF6F"), labels = c("Introgression Tract", "Genome-Wide Background")) +
  labs(x = "", y = "Probability Density", fill = "") +
  theme_blank() +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(title = NULL))

# Print the density plot
print(density_plot)








