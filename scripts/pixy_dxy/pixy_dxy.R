library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/pixy_dxy/")

# 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
introgression_tract_dxy_dist_file_90 = "introgression_tracts_dxy_dist_90.csv"
species_tree_tract_dxy_dist_file_90 = "species_tree_tracts_dxy_dist_90.csv"

dist_introgression_tract_dxy_90 <- scan(introgression_tract_dxy_dist_file_90, what = numeric(), sep = ",")
dist_species_tree_tract_dxy_90 <- scan(species_tree_tract_dxy_dist_file_90, what = numeric(), sep = ",")

num_items = length(dist_introgression_tract_dxy_90)

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df_dxy_90 <- data.frame(
  dist_type = factor(c(rep("introgression_tract", num_items), rep("species_tree_tract", num_items))),
  mean_divergence = c(dist_introgression_tract_dxy_90, dist_species_tree_tract_dxy_90))

mean_values <- aggregate(mean_divergence ~ dist_type, df_dxy_90, mean)

dxy_90 <- ggplot(data = df_dxy_90, aes(x = mean_divergence, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = mean_divergence, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = mean_divergence, label = round(mean_divergence, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgression_tract" = "#1F78B4", "species_tree_tract" = "#333333"),
                     labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")) +
  labs(
    x = "Mean Nucleotide Divergence",
    y = "Probability Density",
    color = "Dist Type"
  ) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) +
  guides(color = guide_legend(title = NULL))


# Display the plot
print(dxy_90)

ggsave(filename = "dxy_90.png", plot = dxy_90)

# 80% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
introgression_tract_dxy_dist_file_80 = "introgression_tracts_dxy_dist_80.csv"
species_tree_tract_dxy_dist_file_80 = "species_tree_tracts_dxy_dist_80.csv"

dist_introgression_tract_dxy_80 <- scan(introgression_tract_dxy_dist_file_80, what = numeric(), sep = ",")
dist_species_tree_tract_dxy_80 <- scan(species_tree_tract_dxy_dist_file_80, what = numeric(), sep = ",")

num_items = length(dist_introgression_tract_dxy_80)

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df_dxy_80 <- data.frame(
  dist_type = factor(c(rep("introgression_tract", num_items), rep("species_tree_tract", num_items))),
  mean_divergence = c(dist_introgression_tract_dxy_80, dist_species_tree_tract_dxy_80))

# Calculate Z-scores for mean_divergence within each dist_type
df_dxy_80_filtered <- df_dxy_80 %>%
  group_by(dist_type) %>%
  mutate(z_score = scale(mean_divergence)) %>%
  filter(abs(z_score) < 4) %>%
  ungroup()

# Remove the z_score column if you no longer need it
df_dxy_80_filtered <- select(df_dxy_80_filtered, -z_score)


mean_values <- aggregate(mean_divergence ~ dist_type, df_dxy_80_filtered, mean)

dxy_80 <- ggplot(data = df_dxy_80_filtered, aes(x = mean_divergence, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = mean_divergence, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = mean_divergence, label = round(mean_divergence, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 3) +
  scale_color_manual(values = c("introgression_tract" = "#1F78B4", "species_tree_tract" = "#333333"),
                     labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")) +
  labs(
    x = "Mean Nucleotide Divergence",
    y = "Probability Density",
    color = "Dist Type"
  ) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) +
  guides(color = guide_legend(title = NULL))


# Display the plot
print(dxy_80)

ggsave(filename = "dxy_80.png", plot = dxy_80)

