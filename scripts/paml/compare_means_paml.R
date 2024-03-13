library(ggplot2)
library(cowplot)

# 90% posterior probability threshold
#-------------------------------------------------------------------------------

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/paml/")

# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dif_dist_file_90 = "difference_distribution_90.csv"
dif_dist_90 <- read.csv(dif_dist_file_90)

paml_dist_file_90 = "paml_dist_90.csv"
paml_dist_90 = read.csv(paml_dist_file_90)

test_stat_dN = mean(paml_dist_90$dN[paml_dist_90$dist_type == "introgressed"]) - mean(paml_dist_90$dN[paml_dist_90$dist_type == "non_introgressed"])

test_stat_dS = mean(paml_dist_90$dS[paml_dist_90$dist_type == "introgressed"]) - mean(paml_dist_90$dS[paml_dist_90$dist_type == "non_introgressed"])

test_stat_dNdS = mean(paml_dist_90$dNdS[paml_dist_90$dist_type == "introgressed"]) - mean(paml_dist_90$dNdS[paml_dist_90$dist_type == "non_introgressed"])
  
# Create a list with response variable names and label expressions
response_vars <- list(
  dN = list(name = "dN", label = expression(italic("d")[N]), test_stat= test_stat_dN),
  dS = list(name = "dS", label = expression(italic("d")[S]), test_stat= test_stat_dS),
  dNdS = list(name = "dNdS", label = expression(italic("d")[N]/italic("d")[S]), test_stat= test_stat_dNdS)
)

# Create empty variables to store the plots
dN_plot <- NULL
dS_plot <- NULL
dNdS_plot <- NULL

# Iterate through response variables
# Iterate through response variables
for (response_var in response_vars) {
  # Create the plot for the current response variable
  current_plot <- ggplot(data = dif_dist_90) +
    stat_density(aes(x = .data[[response_var$name]]), geom="line", position="identity", linewidth = 1.1, color = "#1F78B4") +
    geom_vline(xintercept = response_var$test_stat, linetype = "dashed", color = "darkred", size = 0.7) + 
    labs(title = response_var$label, 
         x = "Difference in Means",
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
  
  # Save the current plot to the corresponding variable
  assign(paste(response_var$name, "_plot", sep = ""), current_plot)
}

# Now you can access the plots using the correct variable names
dN_plot
dS_plot
dNdS_plot


combined_plot <- plot_grid(
  dN_plot + theme(legend.position="none"),
  dS_plot + theme(legend.position="none"),
  dNdS_plot + theme(legend.position="none"),
  labels = c('a.', 'b.', 'c.'),
  nrow=1
  #label_y = 1.05,
  #rel_widths = c(1, 1),
  #rel_heights = c(0.8, 0.80, 0.8)
)

combined_plot
