library("rjson")
library(ggplot2)
library(tidyverse)
library(ggarchery)

setwd("/Users/matt/Documents/Github/dissertation_chapter_2/phylonet_hmm/visualize_tracts/glutathione_peroxidase_1/")

introgression_probabilities = fromJSON(file = "NW_022145611.1.json")
class(introgression_probabilities)

coordinates = read.table("NW_022145611.1_coordinates", header = FALSE)
coordinates = as.vector(coordinates[, 1])

#Original Start: 697300

# Create a data frame with the x and y values
data_df <- data.frame(coordinates = coordinates[697310:697757], 
                      introgression_probabilities = introgression_probabilities[697310:697757])

ggplot(data_df, aes(x = coordinates, y = introgression_probabilities)) +
  geom_line(color = "blue", linetype = "solid") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black") + # Add the dashed line
  xlab("Position on Chromosome NW_022145601.1 (bp)") + ylab("Probability of Introgressed Ancestry") +
  # Remove the background from the theme
  theme(panel.background = element_blank()) +
  # Add lines representing the X and Y Axes
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Rotate x-axis tick labels by 45 degrees
        axis.text.y = element_text(size = 16),  # Keep y-axis tick labels horizontal
        axis.title = element_text(size = 20, face = "bold")) +  # Bold axis titles
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1.00)) +
  # Overlay rectangles and lines representing gene coordinates at the top of the plot
  annotate("rect", xmin = 34800855, xmax = 34801115, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34799748, xmax = 34799818, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34798950, xmax = 34799229, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  geom_arrowsegment(aes(x = 34801230, y = 1.025, xend = 34797264, yend = 1.025), arrow_positions = 0.27, arrows = arrow(length=unit(.35, 'cm'), ends = "last", type = "closed")) + 
  annotate("text", x = mean(c(34797264, 34801230)), y = 1.07, 
           label = "glutathione peroxidase 1", color = "black", vjust = 0) + 
  geom_arrowsegment(aes(x = 34783707, y = 1.025, xend = 34795301, yend = 1.025), arrow_positions = 0.08, arrows = arrow(length=unit(.35, 'cm'), ends = "last", type = "closed")) + 
  annotate("rect", xmin = 34783707, xmax = 34783795, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34785015, xmax = 34785201, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34785686, xmax = 34785769, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34786344, xmax = 34786547, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34787132, xmax = 34787290, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34787820, xmax = 34787978, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34789028, xmax = 34789144, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34789585, xmax = 34789731, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34790724, xmax = 34790921, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34794422, xmax = 34795301, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("text", x = mean(c(34783707, 34795301)), y = 1.07, 
           label = "aarF domain containing kinase 1", color = "black", vjust = 0) + 
  # Add commas representing to the x axis tick marks 
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

ggsave("gpx.svg")
