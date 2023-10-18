library(ggplot2)

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/data/phylonet_hmm/output_visualizations/tracts_by_scaffold/")

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

# Read the scaffold lengths file
scaffold_lengths <- read.table("scaffold_lengths.txt", sep = "\t", header = FALSE)

# Rename the columns
colnames(scaffold_lengths) <- c("chromosome", "length")

# Create a data frame for plotting
plot_data <- data.frame(
  chromosome = factor(genes$chromosome),
  start = genes$start,
  end = genes$end
)

# Create the plot
figure <- ggplot(plot_data) +
  geom_segment(
    aes(x = start, xend = end, y = chromosome, yend = chromosome),
    color = "blue", linewidth = 5
  ) +
  geom_segment(
    data = scaffold_lengths,
    aes(x = 0, xend = length, y = chromosome, yend = chromosome),
    color = "black", size = 0.3
  ) + 
  labs(x = "Position along Chromosome (base pairs)", y = "Chromosome") +
  scale_x_continuous(labels = scales::comma) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 16)
  )

figure

#ggsave(filename = "tracts_80.svg", plot = figure)
ggsave(filename = "tracts_90.svg", plot = figure)
ggsave(filename = "tracts_90.png", plot = figure)
