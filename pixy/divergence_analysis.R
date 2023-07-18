library(sjPlot)
library(fitdistrplus)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)

# Change margin settings
par(mar=c(2,2,2,1))

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/pixy/")

# Define the csv files with list of dXY values 
#introgression_tract_dxy_file = "introgression_tracts_mean_dxy_dist.csv"
introgression_tract_dxy_file = "introgression_tracts_mean_dxy_dist_80.csv"

#species_tree_tract_dxy_file = "species_tree_tracts_mean_dxy.csv"
species_tree_tract_dxy_file = "species_tree_tracts_mean_dxy_dist_80.csv"

# Read the CSV file as a one-dimensional vector
dist_introgression_tract_dxy <- scan(introgression_tract_dxy_file, what = numeric(), sep = ",")
dist_species_tree_tract_dxy <- scan(species_tree_tract_dxy_file, what = numeric(), sep = ",")

# Create data frame
df <- data.frame(
  tract_type = c(rep("introgression_tract", length(dist_introgression_tract_dxy)),
           rep("species_tree_tract", length(dist_species_tree_tract_dxy))),
  divergence = c(dist_introgression_tract_dxy, dist_species_tree_tract_dxy)
)

df_filtered <- df[df$divergence != 0, ]

boxplot(divergence~tract_type,data=df, ylab = "Divergence (Mean dXY)", xlab = "Type")
hist(df$divergence,ylab="Divergence (dXY)", main = "Divergence")

# Plot histogram and overlay normal fit 
fit.norm=fitdist(df$divergence,"norm")
plot(fit.norm)

# Plot histogram and overlay log-normal fit 
fit.lnorm=fitdist(df_filtered$divergence,"lnorm")
plot(fit.lnorm)

# Examine aic for normal fit
fit.norm$aic

# Examine aic for log-normal fit
fit.lnorm$aic

# The data should not be log transformed 

# Check for homogeneity of variances by calculating the variance (sd^2) of each treatment
out=aggregate(divergence~tract_type,data=df,FUN="sd")
out$xvar=out$divergence^2
print(out)

# Run t-test
t.test(divergence~tract_type, data=df, paired=FALSE, var.eq=T)

# Construct bar plot with confidence intervals to show results of t-test

# Use dplyr package to calculate the length, mean, standard deviation, standard error,and confidence intervals of the Density variable for each treatment
summary = df %>%
  group_by(tract_type) %>%
  dplyr::summarize(n=n(),mean=mean(divergence),sd=sd(divergence)) %>%
  mutate(se=sd/sqrt(n)) %>%
  mutate(ic=se * qt((1-0.05)/2 + .5, n-1))

# Examine summary tibble
print(summary)

# Construct plot using summary tibble
ggplot(summary) +
  geom_bar(aes(x=tract_type, y=mean), stat="identity", fill="firebrick") +
  geom_errorbar(aes(x=tract_type, ymin=mean-ic, ymax=mean+ic),
                width=0.4, colour="firebrick4", alpha=0.9, size=1.5)+
  ylab("Mean (CI)")

# The error bar does not overlap the mean of the other treatment in either case. The means are different 

### Density Plot

#Calculate the mean for each group
mu <- ddply(df, "tract_type", summarise, grp.mean=mean(divergence))

figure <- ggplot(df, aes(x=divergence, fill=tract_type)) +
  geom_density(alpha=0.4) +
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=tract_type), show.legend = FALSE, linetype="dashed") + 
  theme(panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_fill_discrete(labels=c('Introgression Tracts', 'Genome-Wide Background')) + theme(legend.title = element_blank())

figure

ggsave(filename = "dXY.png", plot = figure)

ggplot(df, aes(x=divergence, color=tract_type, fill=tract_type)) + 
  geom_histogram(aes(y=..density..), alpha=0.4, 
                 position="identity")+
  geom_density(alpha = 0) +
  theme(panel.background = element_blank()) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_fill_discrete(labels=c('Introgression Tracts', 'Genome-Wide Background')) + theme(legend.title = element_blank())