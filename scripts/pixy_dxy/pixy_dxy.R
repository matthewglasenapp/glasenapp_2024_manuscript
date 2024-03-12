library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/pixy_dxy/")

# 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dxy_dist_file_90 = "dxy_dist_90.csv"
dxy_dist_90 = read.csv(dxy_dist_file_90)

num_items = length(dxy_dist_90$introgression_tract)

df_dxy_90 <- pivot_longer(dxy_dist_90, cols = c(introgression_tract, species_tree_tract), names_to = "dist_type", values_to = "divergence")

df_dxy_90 <- df_dxy_90 %>% arrange(dist_type)

mean_values <- aggregate(divergence ~ dist_type, df_dxy_90, mean)

dxy_90 <- ggplot(data = df_dxy_90, aes(x = divergence, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = divergence, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
 geom_text(data = mean_values, aes(x = divergence, label = round(divergence, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgression_tract" = "#1F78B4", "species_tree_tract" = "#333333"),
                     labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")) +
  labs(
    x = "Absolute Nucleotide Divergence (dXY)",
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
ggsave(filename = "/Users/matt/Desktop/figure_3.pdf", plot = dxy_90)

# 80% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dxy_dist_file_80 = "dxy_dist_80.csv"
dxy_dist_80 = read.csv(dxy_dist_file_80)

num_items = length(dxy_dist_80$introgression_tract)

df_dxy_80 <- pivot_longer(dxy_dist_80, cols = c(introgression_tract, species_tree_tract), names_to = "dist_type", values_to = "divergence")

df_dxy_80 <- df_dxy_80 %>% arrange(dist_type)

df_dxy_80_filtered <- df_dxy_80 %>%
  group_by(dist_type) %>%
  mutate(z_score = scale(divergence)) %>%
  filter(abs(z_score) < 4) %>%
  ungroup() %>%
  select(dist_type, divergence)

mean_values <- aggregate(divergence ~ dist_type, df_dxy_80, mean)

dxy_80 <- ggplot(data = df_dxy_80_filtered, aes(x = divergence, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = divergence, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = divergence, label = round(divergence, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgression_tract" = "#1F78B4", "species_tree_tract" = "#333333"),
                     labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")) +
  labs(
    x = "Absolute Nucleotide Divergence (dXY)",
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

dxy_80
# Display the plot
ggsave(filename = "/Users/matt/Desktop/pixy_80.pdf", plot = dxy_80)

--------
mean(df_dxy_90$divergence[df_dxy_90$dist_type == "introgression_tract"]) 
model1 <- lm(df_dxy_90$divergence[df_dxy_90$dist_type == "species_tree_tract"] ~ 1)
confint(model1, level=0.95)

mean(df_dxy_80$divergence[df_dxy_80$dist_type == "introgression_tract"]) 
model1 <- lm(df_dxy_80$divergence[df_dxy_80$dist_type == "species_tree_tract"] ~ 1)
confint(model1, level=0.95)

