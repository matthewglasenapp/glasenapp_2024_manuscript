library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(ggsignif)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/paml/")

data <- read.csv("paml_dist_80.csv")

# Create a list with response variable names and label expressions
response_vars <- list(
  dN = list(name = "dN", label = expression(italic("d")[N])),
  dS = list(name = "dS", label = expression(italic("d")[S])),
  dNdS = list(name = "dNdS", label = expression(italic("d")[N]/italic("d")[S]))
)

# Create empty variables to store the plots
dN_plot <- NULL
dS_plot <- NULL
dNdS_plot <- NULL

# Iterate through response variables
for (response_var in response_vars) {
  # Create the plot for the current response variable
  current_plot <- ggplot(data, aes(x = dist_type, y = .data[[response_var$name]], fill = dist_type)) +
    geom_boxplot(alpha = 0.3, show.legend = FALSE, width = 0.5, outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = dist_type), alpha = 0.6, show.legend = FALSE, height = 0, width = 0.25, size = 0.5) +
    #geom_signif(comparisons = list(c("introgressed", "non_introgressed")), 
    #map_signif_level=FALSE, annotations = c("")) + 
    #annotate("text", x = c(1.5), y = c(1.5), label = "p = 0.01", size = 3, color = "black") +
    scale_x_discrete(labels = c("introgressed" = "Introgressed", "non_introgressed" = "Non-Introgressed")) +
    labs(x = "", y = response_var$label) +
    scale_fill_manual(values = c("introgressed" = "#1F78B4", "non_introgressed" = "#333333")) +
    scale_color_manual(values = c("introgressed" = "#1F78B4", "non_introgressed" = "#333333")) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.title.y = element_text(size = 10),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
  
  # Save the current plot to the corresponding variable
  assign(paste("plot_", response_var$name, sep = ""), current_plot)
}

# dXY
#-------------------------------------------------------------------------------
setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/pixy_dxy/")

# 80% posterior probability threshold
#-------------------------------------------------------------------------------
# Load the csv files containing the distribution of mean dXY values and read as a one-dimensional vector
dxy_dist_file_80 = "dxy_dist_80.csv"
dxy_dist_80 = read.csv(dxy_dist_file_80)

num_items = length(dxy_dist_80$introgression_tract)

df_dxy_80 <- pivot_longer(dxy_dist_80, cols = c(introgression_tract, species_tree_tract), names_to = "dist_type", values_to = "divergence")

df_dxy_80 <- df_dxy_80 %>% arrange(dist_type)

mean_values <- aggregate(divergence ~ dist_type, df_dxy_80, mean)

dxy_80 <- ggplot(data = df_dxy_80, aes(x = divergence, color = dist_type)) +
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

# Gene Density
#-------------------------------------------------------------------------------
setwd("/Users/matt/Documents/Github/glasenapp_2024_manuscript/data/gene_density/")

# Define the csv files for gene count
gene_count_80_files <- c(
  "introgression_tract_gene_count_80.csv",
  "species_tree_tract_gene_count_80.csv"
)

# Read the CSV files as one-dimensional vectors
dist_gene_counts_80 <- lapply(gene_count_80_files, function(file) {
  scan(file, what = numeric(), sep = ",") / 21.820425
})

# Create data frames for gene count
df_gene <- data.frame(
  gene_count = unlist(dist_gene_counts_80),
  posterior_probability = factor(rep("probability_80", length(unlist(dist_gene_counts_80)))),
  dist_type = factor(rep(c("introgression_tract", "species_tree_tract", "introgression_tract", "species_tree_tract"), each = length(dist_gene_counts_80[[1]])))
)

# Calculate mean and standard deviation for gene counts (only for introgression_tract)
mean_gene <- mean(df_gene[df_gene$dist_type == "introgression_tract",]$gene_count)
sd_gene <- sd(df_gene[df_gene$dist_type == "introgression_tract",]$gene_count)

# Plot gene count data
gene <- ggplot(data = df_gene, aes(x = gene_count)) +
  stat_density(data = df_gene[df_gene$dist_type == "species_tree_tract",], geom="line", position="identity", size = 1.1, color = "#333333") + 
  geom_point(data = data.frame(gene_count = mean_gene, dist_type = "introgression_tract"), aes(x = mean_gene, y = 0.2), color = "#1F78B4", size = 2) +
  geom_errorbarh(data = data.frame(gene_count = mean_gene, dist_type = "introgression_tract"), aes(xmin = mean_gene - sd_gene, xmax = mean_gene + sd_gene, y = 0.2), color = "#1F78B4", height = 0.01, size = 0.5) +
  labs(
    x = "Protein-Coding Genes / Mb",
    y = "Probability Density"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) + 
  guides(color = guide_legend(title = NULL))

# Define the csv files for base count
base_count_80_files <- c(
  "introgression_tract_base_count_80.csv",
  "species_tree_tract_base_count_80.csv"
)

# Read the CSV files as one-dimensional vectors
dist_base_counts_80 <- lapply(base_count_80_files, function(file) {
  scan(file, what = numeric(), sep = ",") / 21820425 * 100
})

# Create data frames for base count
df_base <- data.frame(
  base_count = unlist(dist_base_counts_80),
  posterior_probability = factor(rep("probability_80", length(unlist(dist_base_counts_80)))),
  dist_type = factor(rep(c("introgression_tract", "species_tree_tract", "introgression_tract", "species_tree_tract"), each = length(dist_base_counts_80[[1]])))
)

# Calculate mean and standard deviation for base counts (only for introgression_tract)
mean_base <- mean(df_base[df_base$dist_type == "introgression_tract",]$base_count)
sd_base <- sd(df_base[df_base$dist_type == "introgression_tract",]$base_count)

# Plot base count data
base <- ggplot(data = df_base[df_base$dist_type == "species_tree_tract", ], aes(x = base_count)) +
  stat_density(geom="line", position="identity", size = 1.1, color = "#333333") + 
  geom_point(data = data.frame(base_count = mean_base, dist_type = "introgression_tract"), aes(x = mean_base, y = 1.0), color = "#1F78B4", size = 2) +
  geom_errorbarh(data = data.frame(base_count = mean_base, dist_type = "introgression_tract"), aes(xmin = mean_base - sd_base, xmax = mean_base + sd_base, y = 1.0), color = "#1F78B4", height = 0.08, size = 0.5) +
  labs(
    x = "Percent Coding",
    y = "Probability Density"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) +
  guides(color = guide_legend(title = NULL))

# Define the csv files for psg count
psg_count_80_files <- c(
  "introgression_tract_psg_count_80.csv",
  "species_tree_tract_psg_count_80.csv"
)

# Read the CSV files as one-dimensional vectors
dist_psg_counts_80 <- lapply(psg_count_80_files, function(file) {
  scan(file, what = numeric(), sep = ",") / 21.820425
})

# Create data frames for psg count
df_psg <- data.frame(
  psg_count = unlist(dist_psg_counts_80),
  posterior_probability = factor(rep("probability_80", length(unlist(dist_psg_counts_80)))),
  dist_type = factor(rep(c("introgression_tract", "species_tree_tract"), each = length(dist_psg_counts_80[[1]])))
)

# Calculate mean and standard deviation for psg counts (only for introgression_tract)
mean_psg <- mean(df_psg[df_psg$dist_type == "introgression_tract",]$psg_count)
sd_psg <- sd(df_psg[df_psg$dist_type == "introgression_tract",]$psg_count)

# Plot psg count data
psg <- ggplot(data = df_psg[df_psg$dist_type == "species_tree_tract", ], aes(x = psg_count)) +
  stat_density(geom="line", position="identity", size = 1.1, color = "#333333") + 
  geom_point(data = data.frame(psg_count = mean_psg, dist_type = "introgression_tract"), aes(x = mean_psg, y = 0.85), color = "#1F78B4", size = 2) +
  geom_errorbarh(data = data.frame(psg_count = mean_psg, dist_type = "introgression_tract"), aes(xmin = mean_psg - sd_psg, xmax = mean_psg + sd_psg, y = 0.85), color = "#1F78B4", height = 0.05, size = 0.5) +
  labs(
    x = "Positively Selected Genes / Mb",
    y = "Probability Density"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 10), 
    axis.line = element_line(color = "black", linewidth = 0.5),
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.09))) +
  guides(color = guide_legend(title = NULL))

# Combined plot
#-------------------------------------------------------------------------------
combined_plot <- plot_grid(
  #dxy_80 + theme(legend.position="none"),
  base + theme(legend.position="none"),
  gene + theme(legend.position="none"),
  psg + theme(legend.position="none"),
  plot_dN + theme(legend.position="none"),
  plot_dS + theme(legend.position="none"),
  plot_dNdS + theme(legend.position="none"),
  align = 'vh',
  labels = c('a.', 'b.', 'c.', 'd.', 'e.', 'f.'),
  nrow=2
  #label_y = 1.05,
  #rel_widths = c(1, 1),
  #rel_heights = c(0.8, 0.80, 0.8)
)

combined_plot

ggsave(filename = "/Users/matt/Desktop/nonrandom_80.pdf", plot = combined_plot)

mean(data$dN[data$dist_type == "introgressed"]) 
model1 <- lm(data$dN[data$dist_type == "non_introgressed"] ~ 1)
confint(model1, level=0.95)

mean(data$dS[data$dist_type == "introgressed"]) 
model1 <- lm(data$dS[data$dist_type == "non_introgressed"] ~ 1)
confint(model1, level=0.95)

mean(data$dNdS[data$dist_type == "introgressed"]) 
model1 <- lm(data$dNdS[data$dist_type == "non_introgressed"] ~ 1)
confint(model1, level=0.95)

mean(data$dNdS[data$dist_type == "introgressed"]) 
model1 <- lm(data$dNdS[data$dist_type == "non_introgressed"] ~ 1)
confint(model1, level=0.95)

mean(df_gene$gene_count[df_gene$dist_type == "introgression_tract"])
mean(df_base$base_count[df_base$dist_type == "introgression_tract"])
mean(df_psg$psg_count[df_psg$dist_type == "introgression_tract"])
model5 <- lm(df_gene$gene_count[df_gene$dist_type == "species_tree_tract"] ~ 1)
model6 <- lm(df_base$base_count[df_base$dist_type == "species_tree_tract"] ~ 1)
model7 <- lm(df_psg$psg_count[df_psg$dist_type == "species_tree_tract"] ~ 1)
confint(model5, level=0.95)
confint(model6, level=0.95)
confint(model7, level=0.95)

