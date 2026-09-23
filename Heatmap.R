library(tidyverse)

# 1. Define filename
file_path <- "4562_Summary.csv"

# 2. Load CSV
df <- read_csv(file_path)

# 3. Extract Well ID from Image filename
# Format: {Column}{Row}.tif e.g. "1A.tif", "12H.tif"
df_plot <- df %>%
  select(Image, IntDen_pct_of_max) %>%
  mutate(
    Column = as.numeric(str_extract(Image, "^\\d+")),
    Row    = str_extract(Image, "[A-Z](?=\\.tif)")
  ) %>%
  mutate(Row = factor(Row, levels = rev(sort(unique(Row)))))

# 4. Generate plot title from filename
plot_title <- tools::file_path_sans_ext(basename(file_path))

# 5. Create heatmap
heatmap_plot <- ggplot(df_plot, aes(x = factor(Column), y = Row, fill = IntDen_pct_of_max)) +
  geom_tile(color = "white", lwd = 0.5) +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 100)) +
  coord_fixed() +
  labs(
    title = plot_title,
    x     = "Plate Columns",
    y     = "Plate Rows",
    fill  = "GFP Intensity\n(% of max)"
  ) +
  theme_minimal() +
  theme(
    panel.grid  = element_blank(),
    axis.ticks  = element_blank(),
    plot.title  = element_text(face = "bold", size = 16),
    axis.title  = element_text(face = "italic")
  )

# 6. Save output
ggsave("4562_heatmap_output.png", plot = heatmap_plot, width = 8, height = 6, dpi = 300)

# 7. Bar chart of all wells in descending order with top 12 boxed

# Sort by descending IntDen_pct_of_max and create well label
df_bar <- df_plot %>%
  mutate(Well = paste0(Column, Row)) %>%
  arrange(desc(IntDen_pct_of_max)) %>%
  mutate(
    Rank      = row_number(),
    Top12     = Rank <= 12,
    Well      = factor(Well, levels = Well)  # preserve sort order in plot
  )

# Calculate box boundaries around top 12 bars
box_xmin <- 0.5
box_xmax <- 12.5
box_ymin <- 0
box_ymax <- max(df_bar$IntDen_pct_of_max, na.rm = TRUE) * 1.08

bar_plot <- ggplot(df_bar, aes(x = Well, y = IntDen_pct_of_max, fill = Top12)) +
  geom_col(width = 0.7) +

  # Box around top 12
  annotate("rect",
    xmin  = box_xmin,
    xmax  = box_xmax,
    ymin  = box_ymin,
    ymax  = box_ymax,
    fill  = NA,
    color = "black",
    linewidth = 0.8,
    linetype  = "solid"
  ) +

  # Label the box
  annotate("text",
    x     = (box_xmin + box_xmax) / 2,
    y     = box_ymax,
    label = "Top 12",
    vjust = -0.4,
    size  = 3.5,
    fontface = "bold"
  ) +

  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "steelblue"),
    guide  = "none"
  ) +
  scale_y_continuous(
    limits = c(0, box_ymax * 1.05),
    expand = c(0, 0)
  ) +
  labs(
    title = paste0(plot_title, " — Wells Ranked by GFP Intensity"),
    x     = "Well",
    y     = "GFP Intensity (% of max)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title   = element_text(face = "bold", size = 14),
    axis.title   = element_text(face = "italic")
  )

# Save bar chart
ggsave("GFP_barchart_output.png", plot = bar_plot, width = 14, height = 6, dpi = 300)