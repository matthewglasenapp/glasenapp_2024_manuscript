library(ggplot2)
library(fitdistrplus)

# Change margin settings
par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/count_genes/")

# Define the csv files with list of gene and base countscounts 
introgression_tract_psg_count = "introgression_tract_psg_count.csv"
species_tree_tract_psg_count = "species_tree_tract_psg_count.csv"

# Read the CSV file as a one-dimensional vector
dist_introgression_tract_psg <- scan(introgression_tract_psg_count, what = numeric(), sep = ",")
dist_species_tree_tract_psg <- scan(species_tree_tract_psg_count, what = numeric(), sep = ",")

df <- data.frame(
  tract_type = c(rep("introgression_tract", length(dist_introgression_tract_psg)),
                 rep("species_tree_tract", length(dist_species_tree_tract_psg))),
  psg_count = c(dist_introgression_tract_psg, dist_species_tree_tract_psg))

boxplot(psg_count~tract_type,data=df, ylab = "Overlapping Genes", xlab = "Type")

ggplot(df, aes(x=psg_count, fill=tract_type)) +
  geom_density(alpha=0.4) +
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=tract_type), show.legend = FALSE, linetype="dashed") + 
  theme(panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_fill_discrete(labels=c('Introgression Tracts', 'Genome-Wide Background')) + theme(legend.title = element_blank())

t.test(psg_count~tract_type, data=df, paired=FALSE, var.eq=T)

mean(df$tract_type)
