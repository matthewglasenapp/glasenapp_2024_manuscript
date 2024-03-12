setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/phylonet_hmm/output_visualizations/tract_length_dist/")

library(ggplot2)

dist_90_file = read.csv("tract_dist_90.csv",header=F)
dist_80_file = read.csv("tract_dist_80.csv",header=F)
dist_90 = unlist(dist_90_file)
dist_80 = unlist(dist_80_file)

options(scipen = 999)
figure1 <- ggplot() + 
  aes(dist_90) + 
  geom_histogram(binwidth=10000, colour="black", fill="light blue") + 
  xlab("Length (bases)") + 
  ylab("Number of Introgression Tracts") + 
  theme(axis.title.x = element_text(size=14, face="bold")) + 
  theme(axis.title.y = element_text(size=14, face="bold")) + 
  theme(axis.text.x= element_text(size=12)) + 
  theme(axis.text.y= element_text(size=12)) + 
  theme(axis.line = element_line(colour = "black", size = 1)) + 
  scale_x_continuous(n.breaks = 6, labels = scales::comma, expand = expansion(mult = c(0, .05))) + 
  scale_y_continuous(expand = expansion(mult = c(0, .03))) + 
  theme(panel.background = element_blank()) + 
  stat_bin(binwidth = 10000, aes(y=..count.., label=ifelse(..count..==0,"",..count..)), geom="text", vjust = -.5)

figure1

ggsave(filename = "tract_dist_80.svg", plot = figure)

------------------------------------------------------------------

setwd("/Users/matt/Documents/GitHub/glasenapp_2024_manuscript/data/phylonet_hmm/output_visualizations/tract_length_dist/")
dist_90_file = read.csv("tract_dist_90.csv",header=F)
dist_80_file = read.csv("tract_dist_80.csv",header=F)
dist_90 = unlist(dist_90_file)
dist_80 = unlist(dist_80_file)

options(scipen = 999)  

# Create a combined data frame with reordered levels
combined_data <- data.frame(
  value = c(dist_90, dist_80),
  group = rep(c("dist_90", "dist_80"), times = c(length(dist_90), length(dist_80)))
)

# Modify the levels of the group variable
combined_data$group <- factor(
  combined_data$group,
  levels = c("dist_90", "dist_80"),
  labels = c("Posterior Probability > 90%", "Posterior Probability > 80%")
)

# Create a histogram with facet_wrap
figure2 <- ggplot(data = combined_data, aes(x = value)) +
  geom_histogram(binwidth = 10000, colour = "black", fill = "light blue") +
  labs(x = "Length (bases)", y = "Number of Introgression Tracts") +
  scale_x_continuous(n.breaks = 6, labels = scales::comma, expand = expansion(mult = c(0.005, .05))) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    text = element_text(size = 16),  # Adjust font size for all text elements
    axis.title = element_text(size = 18),  # Adjust font size for axis titles
    strip.text = element_text(size = 18),  # Adjust font size for facet labels
  ) +
  stat_bin(binwidth = 10000, aes(y = ..count.., label = ifelse(..count.. == 0, "", ..count..)), geom = "text", vjust = -0.5) + 
  facet_wrap(~ group, ncol = 1, strip.position = "top", scales = "free_y")  # Stacked facets in one

figure2

ggsave(filename = "/Users/matt/Desktop/tract_length_dist.pdf", plot = figure2, width=169, units = "mm")

#ggsave(filename = "tract_length_distribution.svg", plot = figure2, width=169, units = "mm")
#ggsave(filename = "tract_length_distribution.png", plot = figure2, width=169, height = 200, units = "mm")





