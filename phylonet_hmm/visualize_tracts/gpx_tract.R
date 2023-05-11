library("rjson")
library(ggplot2)

setwd("/Users/matt/desktop/hmm_R_analysis/")

introgression_probabilities = fromJSON(file = "NW_022145611.1.json")
class(introgression_probabilities)

coordinates = read.table("NW_022145611.1_coordinates", header = FALSE)
coordinates = as.vector(coordinates[, 1])

# Create a data frame with the x and y values
data_df <- data.frame(coordinates = coordinates[697507:697757], 
                      introgression_probabilities = introgression_probabilities[697507:697757])

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
  annotate("rect", xmin = 34800855, xmax = 34801115, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34799748, xmax = 34799818, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("rect", xmin = 34798950, xmax = 34799229, ymin = 1.00, ymax = 1.05, 
           fill = "black", color = "black", size = 1) +
  annotate("segment", x = 34797264, y = 1.025, xend = 34801230, yend = 1.025, color = "black", size = 1,
           arrow = arrow(length=unit(.35, 'cm'), ends = "first", type = "open")) +
  annotate("text", x = mean(c(34797264, 34801230)), y = 1.07, 
           label = "glutathione peroxidase 1", color = "black", vjust = 0) + 
  # Add commas representing to the x axis tick marks 
  scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))

ggsave("gpx.svg")