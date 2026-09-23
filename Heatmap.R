library(tidyverse)
library(readxl)

# 1. Define filename (Change this to your actual file path)
file_path <- "4563 Heat Values.xlsx"

# 2. Load data
df <- read_excel(file_path)

# 3. Process Well IDs (Using your exact column names 'Well' and 'Heat Values')
df_plot <- df %>%
  mutate(
    Row = str_extract(Well, "[A-Z]"),
    Column = as.numeric(str_extract(Well, "\\d+"))
  ) %>%
  mutate(Row = factor(Row, levels = rev(sort(unique(Row)))))

# 4. Generate the plot title matching the filename (removes '.xlsx')
plot_title <- tools::file_path_sans_ext(basename(file_path))

# 5. Create Heatmap
heatmap_plot <- ggplot(df_plot, aes(x = factor(Column), y = Row, fill = `Heat Values`)) +
  geom_tile(color = "white", lwd = 0.5) +  
  scale_fill_gradient(low = "blue", high = "red", limits = c(1, 100)) + 
  coord_fixed() + 
  labs(
    title = plot_title,       # Automatically matches the filename
    x = "Plate Columns",      # Custom X-axis title
    y = "Plate Rows",         # Custom Y-axis title
    fill = "Heat Intensity"   # Custom Legend title
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "italic")
  )

# 6. Save the output as a PNG file
# This saves a crisp, high-res image named '4563_heatmap_output.png' in your working folder
ggsave("4563_heatmap_output.png", plot = heatmap_plot, width = 8, height = 6, dpi = 300)