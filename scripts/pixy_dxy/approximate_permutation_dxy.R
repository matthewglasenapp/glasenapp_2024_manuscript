library(ggplot2)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/pixy_dxy/")

# 90% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dif_dist_file_90 = "difference_distribution_90.csv"
dif_dist_90 <- scan(dif_dist_file_90, what = numeric(), sep = ",")

test_stat = -0.00761273

# Create a ggplot object with geom_density directly on the vector
dxy_90 <- ggplot() +
  stat_density(aes(x = dif_dist_90),geom="line", position="identity", linewidth = 1.1, color = "#1F78B4") +
  geom_vline(xintercept = test_stat, linetype = "dashed", color = "darkred", size = 0.7) + 
  labs(x = "Difference in Means",
       y = "Probability Density") + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09)))

dxy_90

# 80% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dif_dist_file_80 = "difference_distribution_80.csv"
dif_dist_80 <- scan(dif_dist_file_80, what = numeric(), sep = ",")

test_stat = -0.007806314842679969

# Create a ggplot object with geom_density directly on the vector
dxy_80 <- ggplot() +
  stat_density(aes(x = dif_dist_80),geom="line", position="identity", linewidth = 1.1, color = "#1F78B4") +
  geom_vline(xintercept = test_stat, linetype = "dashed", color = "darkred", size = 0.7) + 
  labs(x = "Difference in Means",
       y = "Probability Density") + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09)))

dxy_80
