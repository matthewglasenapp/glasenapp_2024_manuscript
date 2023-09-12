library(rjson)
library(ggplot2)
library(tidyverse)
library(ggarchery)

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/scripts/phylonet_hmm/output_visualizations/prominin/")

introgression_probabilities = fromJSON(file = "NW_022145600.1.json")
class(introgression_probabilities)

coordinates = read.table("NW_022145600.1_coordinates", header = FALSE)
coordinates = as.vector(coordinates[, 1])

data_df <- data.frame(coordinates = coordinates[179642:181256], 
                      introgression_probabilities = introgression_probabilities[179642:181256])

figure <- ggplot(data_df, aes(x = coordinates, y = introgression_probabilities)) +
  geom_line(color = "blue", linetype = "solid") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") + # Add the dashed line
  xlab("Position on Chromosome NW_022145600.1 (base pairs)") + ylab("Probability of Introgressed Ancestry") +
  # Remove the background from the theme
  theme(panel.background = element_blank()) +
  # Add lines representing the X and Y Axes
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 20, face = "bold")) +
  coord_cartesian(ylim = c(0, 1.08)) + 
  # Overlay rectangles and lines representing gene coordinates at the top of the plot
  annotate("rect", xmin = 15117808, xmax = 15118014, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) + 
  annotate("rect", xmin = 15119654, xmax = 15119707, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) + 
  annotate("rect", xmin = 15134381, xmax = 	15134730, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5)
  

figure

data_df[data_df$introgression_probabilities > 0.8, ]

