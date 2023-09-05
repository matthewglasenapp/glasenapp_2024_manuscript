library(ggplot2)
library(fitdistrplus)
library(sjPlot)
library(ggpubr)
library(plyr)
library(dplyr)

# Change margin settings
par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/count_genes/")

# Define the csv files with list of gene counts at 90% posterior probability threshold
introgression_tract_base_count_90 = "introgression_tract_base_count_90.csv"
species_tree_tract_base_count_90 = "species_tree_tract_base_count_90.csv"

# Define the csv files with list of gene counts at 80% posterior probability threshold
introgression_tract_base_count_80 = "introgression_tract_base_count_80.csv"
species_tree_tract_base_count_80 = "species_tree_tract_base_count_80.csv"

# Read the CSV file as a one-dimensional vector
dist_introgression_tract_bases_90 <- scan(introgression_tract_base_count_90, what = numeric(), sep = ",")
dist_species_tree_tract_bases_90 <- scan(species_tree_tract_base_count_90, what = numeric(), sep = ",")
dist_introgression_tract_bases_80 <- scan(introgression_tract_base_count_80, what = numeric(), sep = ",")
dist_species_tree_tract_bases_80 <- scan(species_tree_tract_base_count_80, what = numeric(), sep = ",")

total_bases_90_st = 4644060
total_bases_90_it = 4644092
total_bases_80_st = 27343212
total_bases_80_it = 27476503

dist_introgression_tract_bases_90 <- dist_introgression_tract_bases_90 / total_bases_90_it * 100
dist_introgression_tract_bases_80 <- dist_introgression_tract_bases_80 / total_bases_80_it * 100
dist_species_tree_tract_bases_90 <- dist_species_tree_tract_bases_90 / total_bases_90_st * 100
dist_species_tree_tract_bases_80 <- dist_species_tree_tract_bases_80 / total_bases_80_st * 100

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df <- data.frame(
  dist_type = factor(c(rep("introgression_tract", 1000), rep("species_tree_tract", 1000), rep("introgression_tract", 1000), rep("species_tree_tract", 1000))),
  base_count = c(dist_introgression_tract_bases_90, dist_species_tree_tract_bases_90, dist_introgression_tract_bases_80, dist_species_tree_tract_bases_80),
  posterior_probability = factor(c(rep("probability_90", 2000), rep("probability_80", 2000)))
)

# Reorder levels of posterior_probability
df$posterior_probability <- factor(df$posterior_probability, levels = c("probability_90", "probability_80"))

# Plot data
fig1 <- ggplot(data = df, aes(x = base_count, fill = dist_type)) +
  geom_density(aes(y=after_stat(density))) + 
  facet_wrap(
    ~ posterior_probability,
    ncol = 1,
    scales = "free",
    #scales = "free_y",
    #scales = "free_x",
    labeller = labeller(
      posterior_probability = c(
        "probability_90" = "Posterior Probability Threshold: 90%",
        "probability_80" = "Posterior Probability Threshold: 80%"
      )
    )
  ) +
  scale_fill_manual(values = c("introgression_tract" = "#A6CEE3", "species_tree_tract" = "#FDBF6F"), labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")
  ) + labs(
    x = "Mean Percent Coding",
    y = "Probability Density",
    fill = "Dist Type"
  ) + theme_blank() +
  theme(
    axis.title = element_text(size = 14),
    legend.position = "right",
    panel.grid = element_blank(), 
  ) + 
  guides(fill = guide_legend(title = NULL))

fig1

ggsave(filename = "dist_coding.eps", plot = fig1)
ggsave(filename = "dist_coding.png", plot = fig1)

###

model <- lm(dist_species_tree_tract_bases_80 ~ 1)

# Calculate the confidence interval
confint(model, level=0.95)
