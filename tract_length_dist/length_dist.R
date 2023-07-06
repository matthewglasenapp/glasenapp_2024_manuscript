setwd("/Users/matt/Documents/Github/dissertation_chapter_2/tract_length_dist/")

library(ggExtra)
library(ggplot2)

#x = read.csv("tract_dist.csv",header=F)
#x = read.csv("all_tracts.csv",header=F)
x = read.csv("tract_dist_80.csv",header=F)
y=unlist(x)

options(scipen = 999)
figure <- ggplot() + 
  aes(y) + 
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

figure

ggsave(filename = "tract_dist_80.svg", plot = figure)




