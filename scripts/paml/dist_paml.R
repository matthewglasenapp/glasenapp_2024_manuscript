library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/paml/")

# dN 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dN_introgressed_90_file = "dN_introgressed_90.csv"
dN_non_introgressed_90_file = "dN_non_introgressed_90.csv"

dist_dN_introgressed_90 <- scan(dN_introgressed_90_file, what = numeric(), sep = ",")
dist_dN_non_introgressed_90 <- scan(dN_non_introgressed_90_file, what = numeric(), sep = ",")

num_items = length(dist_dN_introgressed_90)

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df_dN_90 <- data.frame(
  dist_type = factor(c(rep("introgressed", num_items), rep("non_introgressed", num_items))),
  dN = c(dist_dN_introgressed_90, dist_dN_non_introgressed_90))

# Calculate Z-scores for mean_divergence within each dist_type
df_dN_90_filtered <- df_dN_90 %>%
  group_by(dist_type) %>%
  mutate(z_score = scale(dN)) %>%
  filter(abs(z_score) < 4) %>%
  ungroup()

# Remove the z_score column if you no longer need it
df_dN_90_filtered <- select(df_dN_90_filtered, -z_score)

mean_values <- aggregate(dN ~ dist_type, df_dN_90_filtered, mean)

dN_90 <- ggplot(data = df_dN_90_filtered, aes(x = dN, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = dN, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = dN, label = round(dN, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgressed" = "#1F78B4", "non_introgressed" = "#333333"),
                     labels = c("introgressed" = "Introgressed Genes", "non_introgressed" = "Non-Introgressed Genes")) +
  labs(
    x = "dN",
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
print(dN_90)

# dS 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dS_introgressed_90_file = "dS_introgressed_90.csv"
dS_non_introgressed_90_file = "dS_non_introgressed_90.csv"

dist_dS_introgressed_90 <- scan(dS_introgressed_90_file, what = numeric(), sep = ",")
dist_dS_non_introgressed_90 <- scan(dS_non_introgressed_90_file, what = numeric(), sep = ",")

num_items = length(dist_dS_introgressed_90)

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df_dS_90 <- data.frame(
  dist_type = factor(c(rep("introgressed", num_items), rep("non_introgressed", num_items))),
  dS = c(dist_dS_introgressed_90, dist_dS_non_introgressed_90))

# Calculate Z-scores for mean_divergence within each dist_type
df_dS_90_filtered <- df_dS_90 %>%
  group_by(dist_type) %>%
  mutate(z_score = scale(dS)) %>%
  filter(abs(z_score) < 4) %>%
  ungroup()

# Remove the z_score column if you no longer need it
df_dS_90_filtered <- select(df_dS_90_filtered, -z_score)

mean_values <- aggregate(dS ~ dist_type, df_dS_90_filtered, mean)

dS_90 <- ggplot(data = df_dS_90_filtered, aes(x = dS, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = dS, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = dS, label = round(dS, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgressed" = "#1F78B4", "non_introgressed" = "#333333"),
                     labels = c("introgressed" = "Introgressed Genes", "non_introgressed" = "Non-Introgressed Genes")) +
  labs(
    x = "dS",
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
print(dS_90)

# dNdS 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dNdS_introgressed_90_file = "dNdS_introgressed_90.csv"
dNdS_non_introgressed_90_file = "dNdS_non_introgressed_90.csv"

dist_dNdS_introgressed_90 <- scan(dNdS_introgressed_90_file, what = numeric(), sep = ",")
dist_dNdS_non_introgressed_90 <- scan(dNdS_non_introgressed_90_file, what = numeric(), sep = ",")

num_items = length(dist_dNdS_introgressed_90)

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df_dNdS_90 <- data.frame(
  dist_type = factor(c(rep("introgressed", num_items), rep("non_introgressed", num_items))),
  dNdS = c(dist_dNdS_introgressed_90, dist_dNdS_non_introgressed_90))

# Calculate Z-scores for mean_divergence within each dist_type
df_dNdS_90_filtered <- df_dNdS_90 %>%
  group_by(dist_type) %>%
  mutate(z_score = scale(dNdS)) %>%
  filter(abs(z_score) < 4) %>%
  ungroup()

# Remove the z_score column if you no longer need it
df_dNdS_90_filtered <- select(df_dNdS_90_filtered, -z_score)

mean_values <- aggregate(dNdS ~ dist_type, df_dNdS_90_filtered, mean)

dNdS_90 <- ggplot(data = df_dNdS_90_filtered, aes(x = dNdS, color = dist_type)) +
  stat_density(geom="line",position="identity", size = 1.1) + 
  geom_vline(data = mean_values, aes(xintercept = dNdS, linetype = dist_type, color = dist_type),
             linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_text(data = mean_values, aes(x = dNdS, label = round(dNdS, 3), y = 0.0, color = dist_type), vjust = 0.0, hjust = -0.1, show.legend = FALSE, size = 2.5) +
  scale_color_manual(values = c("introgressed" = "#1F78B4", "non_introgressed" = "#333333"),
                     labels = c("introgressed" = "Introgressed Genes", "non_introgressed" = "Non-Introgressed Genes")) +
  labs(
    x = "dNdS",
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
print(dNdS_90)

plot <- plot_grid(
  dN_90 + theme(legend.position="none"),
  dS_90 + theme(legend.position="none"),
  dNdS_90 + theme(legend.position="none"),
  align = 'vh',
  labels = c('a.', 'b.', 'c.'),
  nrow = 1
)

plot

