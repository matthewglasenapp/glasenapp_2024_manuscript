library(rjson)
library(ggplot2)
library(tidyverse)
library(ggarchery)

setwd("/Users/matt/Documents/Github/glasenapp_2024_manuscript/data/phylonet_hmm/output_visualizations/arachidonate_5_lipoxygenase/")

introgression_probabilities = fromJSON(file = "NW_022145601.1.json")
class(introgression_probabilities)

coordinates = read.table("NW_022145601.1_coordinates.txt", header = FALSE)
coordinates = as.vector(coordinates[, 1])

# Create a data frame with the x and y values
data_df <- data.frame(coordinates = coordinates[100:1200], 
                      introgression_probabilities = introgression_probabilities[100:1200])

figure <- ggplot(data_df, aes(x = coordinates, y = introgression_probabilities)) +
  geom_line(color = "blue", linetype = "solid") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black") + # Add the dashed line
  xlab("Position on Chromosome NW_022145601.1 (base pairs)") + ylab("Probability of Introgressed Ancestry") +
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
  annotate("rect", xmin = 170756, xmax = 172402, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 175704, xmax = 175748, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 176657, xmax = 176785, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 177249, xmax = 177346, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 177643, xmax = 177775, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 178403, xmax = 178495, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 179302, xmax = 179516, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 179896, xmax = 180155, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 180545, xmax = 180645, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 181381, xmax = 181554, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 182143, xmax = 182306, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 183906, xmax = 184101, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  geom_arrowsegment(aes(x = 170401, y = 1.025, xend = 184101, yend = 1.025), arrow_positions = 0.3, arrows = arrow(length=unit(.35, 'cm'), ends = "last", type = "closed")) + 
  annotate("text", x = mean(c(170401, 184101)), y = 1.06, 
           label = "arachidonate\n5-lipoxygenase", color = "black", vjust = 0) + 
  geom_arrowsegment(aes(x = 151263, y = 1.025, xend = 162808, yend = 1.025), arrow_positions = 0.3, arrows = arrow(length=unit(.35, 'cm'), ends = "last", type = "closed")) + 
  annotate("text", x = mean(c(151263, 162808)), y = 1.06, 
           label = "lysosomal amino acid\ntransporter 1", color = "black", vjust = 0) + 
  annotate("rect", xmin = 151263, xmax = 151397, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 152335, xmax = 152661, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 153562, xmax = 153688, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 155242, xmax = 155350, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 155931, xmax = 156133, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 156698, xmax = 156784, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 157672, xmax = 157814, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 159949, xmax = 159976, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 160821, xmax = 161039, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 161457, xmax = 161648, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 162063, xmax = 162808, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  geom_arrowsegment(aes(x = 198738, y = 1.025, xend = 210294, yend = 1.025), arrow_positions = 0.33, arrows = arrow(length=unit(.35, 'cm'), ends = "last", type = "closed")) + 
  annotate("text", x = mean(c(198738, 210294)), y = 1.06, 
           label = "glycosyltransferase 8\ndomain-containing protein 1", color = "black", vjust = 0) +
  annotate("rect", xmin = 198738, xmax = 199186, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 201299, xmax = 201391, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 202975, xmax = 203176, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 203916, xmax = 204033, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 204438, xmax = 204626, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 205330, xmax = 205496, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 205792, xmax = 205904, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  annotate("rect", xmin = 206637, xmax = 210294, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 0.5) +
  # Add commas representing to the x axis tick marks 
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) + 
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1.0))

figure

#ggsave(filename="/Users/matt/Desktop/figure_5.pdf",plot=figure, width = 11, height = 8.5)
  
#ggsave(filename = "arach.svg", width = 11, height = 9, plot = figure, dpi = 700)
#ggsave(filename = "arach.png", width = 11, height = 8.5, plot = figure)

