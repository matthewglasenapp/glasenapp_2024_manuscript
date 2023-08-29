library(ggplot2)
library(fitdistrplus)

# Change margin settings
par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/count_genes/")

# Define the csv files with list of gene and base countscounts 
introgression_tract_gene_count = "introgression_tract_gene_count.csv"
introgression_tract_base_count = "introgression_tract_base_count.csv"

species_tree_tract_gene_count = "species_tree_tract_gene_count.csv"
species_tree_tract_base_count = "species_tree_tract_base_count.csv"

# Read the CSV file as a one-dimensional vector
dist_introgression_tract_genes <- scan(introgression_tract_gene_count, what = numeric(), sep = ",")
dist_species_tree_tract_genes <- scan(species_tree_tract_gene_count, what = numeric(), sep = ",")
dist_introgression_tract_bases <- scan(introgression_tract_base_count, what = numeric(), sep = ",")
dist_species_tree_tract_bases <- scan(species_tree_tract_base_count, what = numeric(), sep = ",")

df <- data.frame(
  tract_type = c(rep("introgression_tract", length(dist_introgression_tract_genes)),
                 rep("species_tree_tract", length(dist_species_tree_tract_genes))),
  gene_count = c(dist_introgression_tract_genes, dist_species_tree_tract_genes),
  coding_bases = c(dist_introgression_tract_bases, dist_species_tree_tract_bases)
)

boxplot(gene_count~tract_type,data=df, ylab = "Overlapping Genes", xlab = "Type")
boxplot(coding_bases~tract_type,data=df, ylab = "Overlapping Genes", xlab = "Type")

hist(df$gene_count,ylab="Overlapping Genes", main = "")
hist(df$coding_bases,ylab="Coding Bases", main = "")

# Plot histogram and overlay normal fit 
fit.norm_genes=fitdist(df$gene_count,"norm")
plot(fit.norm_genes)

fit.norm_bases=fitdist(df$coding_bases,"norm")
plot(fit.norm_bases)

# Plot histogram and overlay log-normal fit 
fit.lnorm_genes=fitdist(df$gene_count,"lnorm")
plot(fit.lnorm_genes)

fit.lnorm_bases=fitdist(df$coding_bases,"lnorm")
plot(fit.lnorm_bases)

# Examine aic for normal fit
fit.norm_genes$aic
fit.norm_bases$aic

# Examine aic for log-normal fit
fit.lnorm_genes$aic
fit.lnorm_bases$aic

# The data should not be log transformed 

# Check for homogeneity of variances by calculating the variance (sd^2) of each treatment
out_genes=aggregate(gene_count~tract_type,data=df,FUN="sd")
out_genes$xvar=out_genes$gene_count^2
print(out_genes)

out_bases=aggregate(coding_bases~tract_type,data=df,FUN="sd")
out_bases$xvar=out_bases$coding_bases^2
print(out_bases)


ggplot(df, aes(x=gene_count, fill=tract_type)) +
  geom_density(aes(y=after_stat(density)), alpha=0.4) +
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=tract_type), show.legend = FALSE, linetype="dashed") + 
  theme(panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_fill_discrete(labels=c('Introgression Tracts', 'Genome-Wide Background')) + theme(legend.title = element_blank())

ggplot(df, aes(x=coding_bases, fill=tract_type)) +
  geom_density(alpha=0.4) +
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=tract_type), show.legend = FALSE, linetype="dashed") + 
  theme(panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_fill_discrete(labels=c('Introgression Tracts', 'Genome-Wide Background')) + theme(legend.title = element_blank())

t.test(gene_count~tract_type, data=df, paired=FALSE, var.eq=T)
t.test(coding_bases~tract_type, data=df, paired=FALSE, var.eq=T)



