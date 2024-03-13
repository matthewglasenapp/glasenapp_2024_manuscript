library(ggplot2)
library(viridisLite)
library(svglite)
library(dplyr)

setwd("/Users/matt/Documents/Github/glasenapp_2024_manuscript/data/phylonet_hmm/output_visualizations/tracts_by_scaffold/")

# Read data from CSV file
gene_density_data <- read.csv("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/gene_density/chromosome_window_density.csv", header = FALSE, col.names = c("chromosome", "start", "end", "gene_density"))

# Create a data frame for plotting
df_gene_density <- gene_density_data %>%
  mutate(chromosome = factor(chromosome, levels = unique(chromosome)))

# Calculate window centers
df_gene_density$window_center <- (df_gene_density$start + df_gene_density$end) / 2

mytheme <- theme_classic() + theme(
  legend.position = "top",
  axis.line = element_line(color = "black", size = 1.5),
  axis.text = element_text(color = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1),
  text = element_text(size = 10, face = "bold", color = "black")
)
theme_set(mytheme)

# Read the bed file
#bed_file <- "ten_kb_tracts_80.bed"
bed_file <- "ten_kb_tracts.bed"
genes <- read.table(bed_file, sep = "\t", header = FALSE)
# Rename the columns
colnames(genes) <- c("chromosome", "start", "end", "gene_name")

# Read the gaps file
gap_file <- "100kb_gaps.bed"
gaps <- read.table(gap_file, sep = "\t", header = FALSE)
colnames(gaps) <- c("chromosome", "start", "end", "name")

# Read the scaffold lengths file
scaffold_lengths <- read.table("scaffold_info.txt", sep = "\t", header = FALSE)
# Rename the columns
colnames(scaffold_lengths) <- c("chromosome", "length", "gene_density")
# Order scaffold_lengths by gene_density in descending order
scaffold_lengths <- scaffold_lengths[order(scaffold_lengths$gene_density),]
scaffold_lengths$chromosome <- factor(scaffold_lengths$chromosome, levels = scaffold_lengths$chromosome)

# Create a new column in genes dataframe to represent the order of chromosomes
genes$order <- match(genes$chromosome, scaffold_lengths$chromosome)
# Sort genes by the order of chromosomes
genes <- genes[order(genes$order),]
# Remove the 'order' column if no longer needed
genes$order <- NULL


# Create a data frame for plotting
plot_data <- data.frame(
  chromosome = factor(genes$chromosome, levels = scaffold_lengths$chromosome),
  start = genes$start,
  end = genes$end
)

# Create the plot
figure <- ggplot(plot_data) +
  geom_blank(data = scaffold_lengths, aes(x = 0, y = chromosome)) + 
  geom_tile(data = df_gene_density, aes(x = window_center, y = chromosome, fill = gene_density),
            width = df_gene_density$end - df_gene_density$start + 1, height = 0.2) +
  geom_segment(
    aes(x = start, xend = end, y = chromosome, yend = chromosome),
    color = "black", linewidth = 5
  ) +
  labs(x = "Position along Scaffold (base pairs)", y = "Scaffold") +
  scale_x_continuous(labels = scales::comma) +
  scale_fill_viridis_c(name = "Gene Density", option = "rocket") +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 16),
    legend.position = "inside",
    legend.justification = "right",
    #legend.key.size = unit(0.5, "cm"),
  )

figure

#ggsave(filename = "/Users/matt/Desktop/figure_2.pdf", plot = figure)

#ggsave(filename = "tracts_80.svg", plot = figure)
#ggsave(filename = "/Users/matt/Desktop/tracts_80.svg", plot = figure)
#ggsave(filename = "/Users/matt/Desktop/tracts_90.png", plot = figure)