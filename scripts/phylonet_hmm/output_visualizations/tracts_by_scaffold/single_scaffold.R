library(ggplot2)

setwd("/Users/matt/Documents/GitHub/dissertation_chapter_2/phylonet_hmm/visualize_tracts/tracts_by_scaffold/")

mytheme <- theme_classic() + theme(
  legend.position = "top",
  axis.line = element_line(color = "black", size = 1.5),
  axis.text = element_text(color = "black"),
  axis.text.x = element_text(angle = 45, hjust = 1),
  text = element_text(size = 10, face = "bold", color = "black")
)
theme_set(mytheme)

# Read the bed file
bed_file <- "ten_kb_tracts_80.bed"
genes <- read.table(bed_file, sep = "\t", header = FALSE)

# Rename the columns
colnames(genes) <- c("chromosome", "start", "end", "gene_name")

# Filter the data for the desired chromosome
chromosome_to_plot <- "NW_022145601.1"
plot_data <- subset(genes, chromosome == chromosome_to_plot)

# Read the scaffold lengths file
scaffold_lengths <- read.table("scaffold_lengths.txt", sep = "\t", header = FALSE)

# Rename the columns
colnames(scaffold_lengths) <- c("chromosome", "length")

# Get the length of the chromosome from scaffold lengths data
chromosome_length <- subset(scaffold_lengths, chromosome == chromosome_to_plot)$length

# Create the plot
figure <- ggplot() +
  geom_segment(
    data = plot_data,
    aes(x = start, xend = end, y = 0, yend = 0),
    color = "blue", size = 10
  ) +
  geom_segment(
    data = subset(scaffold_lengths, chromosome == chromosome_to_plot),
    aes(x = 0, xend = chromosome_length, y = 0, yend = 0),
    color = "black", size = 0.5
  ) +
  labs(x = "Position") +
  scale_x_continuous(labels = scales::comma) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 20),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )

figure

# Create the plot
figure <- ggplot() +
  geom_segment(
    data = subset(plot_data, start >= 150000 & end <= 205000),
    aes(x = start, xend = end, y = 0, yend = 0),
    color = "blue", size = 10
  ) +
  geom_segment(
    data = subset(scaffold_lengths, chromosome == chromosome_to_plot),
    aes(x = 0, xend = chromosome_length, y = 0, yend = 0),
    color = "black", size = 0.5
  ) +
  labs(x = "Position") +
  scale_x_continuous(labels = scales::comma, limits = c(150000, 205000)) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 1),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 20),
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.y = element_blank(),
    axis.line.y = element_blank()
  )

figure


# Print the plot
print(figure)


# Print the plot
ggsave(filename = "single_tract_80.svg", plot = figure)

