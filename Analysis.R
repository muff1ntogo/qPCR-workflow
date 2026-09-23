# Install packages if you haven't already:
# install.packages("tidyverse")

library(dplyr)
library(readr)

# Load the data
df <- read_csv("4566 Results.csv")

# Group and aggregate
summary_df <- df %>%
  group_by(Image) %>%
  summarise(
    Colony_Count  = sum(!is.na(IntDen)), # Counts non-missing values (matches pandas .count())
    Total_Area    = sum(Area_px2, na.rm = TRUE),
    Avg_Mean_GFP  = mean(Mean_GFP, na.rm = TRUE),
    Avg_IntDen    = mean(IntDen, na.rm = TRUE),
    Avg_RawIntDen = mean(RawIntDen, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Normalise Avg_IntDen as % of the maximum
    IntDen_pct_of_max = round((Avg_IntDen / max(Avg_IntDen, na.rm = TRUE)) * 100, 2)
  )

# Save to CSV (equivalent to index=False)
write_csv(summary_df, "4566_Summary.csv")

# Print summary to console
print(summary_df)