# Install and load the necessary package
library(ggplot2)

# Create a data frame with the given data
data <- data.frame(
  Category = c("Coding", "Intron", "Intergenic"),
  Bases = c(842951, 3455534, 4773694),
  Percentage = c(9.3, 38.1, 52.6)
)

# Create a pie chart using ggplot2
ggplot(data, aes(x = "", y = Percentage, fill = Category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Category, "\n ", Percentage, "%")), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Blues") +
  guides(fill = FALSE) +
  theme_void()

