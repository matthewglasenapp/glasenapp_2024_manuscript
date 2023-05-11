library("rjson")
library(ggplot2)
library(tidyverse)
library(ggarchery)

setwd("/Users/matt/desktop/hmm_R_analysis/")

introgression_probabilities = fromJSON(file = "NW_022145601.1.json")
class(introgression_probabilities)

coordinates = read.table("NW_022145601.1_coordinates", header = FALSE)
coordinates = as.vector(coordinates[, 1])

# Create a data frame with the x and y values
data_df <- data.frame(coordinates = coordinates[323:930], 
                      introgression_probabilities = introgression_probabilities[323:930])

ggplot(data_df, aes(x = coordinates, y = introgression_probabilities)) +
  geom_line(color = "red", linetype = "solid") +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black") + # Add the dashed line
  xlab("Position on Scaffold NW_022145601.1") + ylab("Probability of Hybrid Ancestry") +
  # Remove the background from the theme
  theme(panel.background = element_blank()) +
  # Add lines representing the X and Y Axes
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", size = 1),
        axis.text = element_text(size = 16), # increase font size of axis labels
        axis.title = element_text(size = 20)) + # increase font size of axis titles
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
  annotate("text", x = mean(c(170401, 184101)), y = 1.07, 
           label = "arachidonate 5-lipoxygenase", color = "black", vjust = 0) + 
  # Add commas representing to the x axis tick marks 
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

ggsave("arach.svg")
