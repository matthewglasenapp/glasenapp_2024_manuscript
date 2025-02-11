library(ggplot2)
library(scales)

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/introgression_by_scaffold_density/")

file = "introgression_by_scaffold_updated_90.tsv"
#file = "introgression_by_scaffold_updated_80.tsv"

data = read.csv(file, sep = "\t")

data$percent_introgressed <- data$combined_tract_length / data$length.spanned.by.first.and.last.sites..bp.*100

data$length_mb <- data$length..bp. / 1000000

fig1 <- ggplot(data, aes(x=length_mb, y=percent_introgressed)) + 
  geom_point(size = 3) +
  labs(x="Chromosome Length (Mb)", y = "Percent Introgressed") + 
  theme_classic() + 
  geom_smooth(method="lm", color = "black") +
  theme(axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold"), axis.text.x = element_text(size=18), axis.text.y = element_text(size=18))

fig1

lm <- lm(percent_introgressed ~ length.spanned.by.first.and.last.sites..bp., data = data)

summary(lm)

#ggsave(filename = "introgression_by_length.svg", plot = fig1, width=169, units = "mm")
#ggsave(filename = "introgression_by_length.png", plot = fig1, width=169, units = "mm")

#-------------
# Gene Density

data$percent_introgressed <- data$combined_tract_length_10kb / data$length.spanned.by.first.and.last.sites..bp.*100

fig2 <- ggplot(data, aes(x=gene_density, y=percent_introgressed)) + 
  geom_point(size = 3) +
  labs(x="Gene Density", y = "Percent Introgressed") + 
  theme_classic() + 
  geom_smooth(method="lm", color = "black") +
  theme(axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold"), axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) + 
  scale_y_continuous(limits = c(0, NA), oob = scales::oob_keep)

fig2

lm <- lm(percent_introgressed ~ gene_density, data = data)

summary(lm)

ggsave(filename = "/Users/matt/Desktop/figure_s3.pdf", plot = fig2, width=169, units = "mm")
#ggsave(filename = "introgression_by_gene_density.png", plot = fig2, width=169, units = "mm")
