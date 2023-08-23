library(ggplot2)

setwd("/Users/matt/Documents/GitHub/dissertation_chapter_2/introgression_by_scaffold_length/")

file = "introgression_by_scaffold_80.csv"

data = read.csv(file)
data$percent_sites_introgressed <- data$number.posterior.probabilities...threshold / data$number.nexus.sites.in.phmm.alignment

ggplot(data, aes(x=length.spanned.by.first.and.last.sites..bp., y=percent_sites_introgressed)) + 
  geom_point(size = 3) +
  labs(x="Chromosome Length", y = "Percent Introgressed") + 
  theme_classic() + 
  geom_smooth(method="lm", color = "black") +
  theme(axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold"), axis.text.x = element_text(size=18), axis.text.y = element_text(size=18))

lm <- lm(percent_sites_introgressed ~ length.spanned.by.first.and.last.sites..bp., data = data)

summary(lm)
