library(ggplot2)

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/data/gene_density/")

# Define the csv files with list of gene counts at 90% posterior probability threshold
introgression_tract_gene_count_90 = "introgression_tract_gene_count_90.csv"
species_tree_tract_gene_count_90 = "species_tree_tract_gene_count_90.csv"

# Define the csv files with list of gene counts at 80% posterior probability threshold
introgression_tract_gene_count_80 = "introgression_tract_gene_count_80.csv"
species_tree_tract_gene_count_80 = "species_tree_tract_gene_count_80.csv"

# Read the CSV file as a one-dimensional vector
dist_introgression_tract_genes_90 <- scan(introgression_tract_gene_count_90, what = numeric(), sep = ",")
dist_species_tree_tract_genes_90 <- scan(species_tree_tract_gene_count_90, what = numeric(), sep = ",")
dist_introgression_tract_genes_80 <- scan(introgression_tract_gene_count_80, what = numeric(), sep = ",")
dist_species_tree_tract_genes_80 <- scan(species_tree_tract_gene_count_80, what = numeric(), sep = ",")

mb_introgressed_90 = 3.747392
mb_introgressed_80 = 21.820425
mb_species_tree_90 = 3.747392
mb_species_tree_80 = 21.820425

dist_introgression_tract_genes_90 <- dist_introgression_tract_genes_90 / mb_introgressed_90
dist_species_tree_tract_genes_90 <- dist_species_tree_tract_genes_90 / mb_species_tree_90
dist_introgression_tract_genes_80 <- dist_introgression_tract_genes_80 / mb_introgressed_80
dist_species_tree_tract_genes_80 <- dist_species_tree_tract_genes_80 / mb_species_tree_80

# Add one-dimensional vectors of mean dXY values to a dataframe called df
df <- data.frame(
  dist_type = factor(c(rep("introgression_tract", 1000), rep("species_tree_tract", 1000), rep("introgression_tract", 1000), rep("species_tree_tract", 1000))),
  gene_count = c(dist_introgression_tract_genes_90, dist_species_tree_tract_genes_90, dist_introgression_tract_genes_80, dist_species_tree_tract_genes_80),
  posterior_probability = factor(c(rep("probability_90", 2000), rep("probability_80", 2000)))
)

# Reorder levels of posterior_probability
df$posterior_probability <- factor(df$posterior_probability, levels = c("probability_90", "probability_80"))

# Plot data
fig1 <- ggplot(data = df, aes(x = gene_count, fill = dist_type)) +
  geom_density(aes(y=after_stat(density)),alpha=0.8) + 
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
  scale_fill_manual(values = c("introgression_tract" = "#A6CEE3", "species_tree_tract" = "#999999"), labels = c("introgression_tract" = "Introgression Tracts", "species_tree_tract" = "Genome-Wide Background")
  ) + labs(
    x = "Protein-Coding Genes / Mb",
    y = "Probability Density",
    fill = "Dist Type"
  ) + theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 12),
    panel.grid = element_blank(), 
    legend.position = "top",
    axis.text = element_text(size = 12)
  ) + 
  guides(fill = guide_legend(title = NULL))

fig1

ggsave(filename = "dist_genes.eps", plot = fig1)
ggsave(filename = "dist_genes.svg", plot = fig1)
ggsave(filename = "dist_genes.png", width=169, units = "mm", plot = fig1)

######

model1 <- lm(dist_species_tree_tract_genes_80 ~ 1)
model2 <- lm(dist_species_tree_tract_genes_90 ~ 1)

# Calculate the confidence interval
confint(model1, level=0.95)
confint(model2, level=0.95)

mean(dist_introgression_tract_genes_80)
mean(dist_introgression_tract_genes_90)

