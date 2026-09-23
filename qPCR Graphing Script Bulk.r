# Install required packages if not already installed
# install.packages(c("readr", "dplyr", "ggplot2", "openxlsx", "tidyr", "scales"))

library(readr)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)
library(scales)

REFERENCE_GENE <- "Gapdh"

# ---- Custom y-axis transform ----
make_squish_trans <- function(squish_break = 10, squish_factor = 6) {
  trans <- function(x) ifelse(x <= squish_break, x, squish_break + (x - squish_break) / squish_factor)
  inv   <- function(x) ifelse(x <= squish_break, x, squish_break + (x - squish_break) * squish_factor)
  scales::trans_new(
    name = paste0("squish_above_", squish_break),
    transform = trans,
    inverse = inv
  )
}

# ---- Read + clean a single qPCR CSV file ----
read_qpcr_csv <- function(file_path) {
  df <- read_csv(file_path, col_types = cols(.default = "c"), show_col_types = FALSE)
  names(df)[1] <- "index"
  df <- df %>% mutate(across(everything(), trimws))

  df$Dox <- dplyr::case_when(
    grepl("\\+", df$Dox) | toupper(df$Dox) == "B" ~ "B",
    grepl("-", df$Dox)   | toupper(df$Dox) == "A" ~ "A",
    TRUE ~ df$Dox
  )

  df$Ct <- suppressWarnings(as.numeric(df$Ct))
  df <- df %>% filter(!is.na(Cell_line), Cell_line != "", !is.na(Primer), Primer != "")
  df$batch <- tools::file_path_sans_ext(basename(file_path))
  df
}

# ---- Calculate ddCt-derived two_fold and Dox+/Dox- ratio ----
calculate_ddct <- function(all_data, reference_gene = REFERENCE_GENE) {
  ref_ct <- all_data %>%
    filter(Primer == reference_gene) %>%
    select(batch, Cell_line, subnumber, Dox, Ct_ref = Ct)

  targets <- all_data %>%
    filter(Primer != reference_gene) %>%
    left_join(ref_ct, by = c("batch", "Cell_line", "subnumber", "Dox")) %>%
    mutate(
      dCt = Ct - Ct_ref,
      two_fold = 2 ^ (-dCt)
    )

  dox_a <- targets %>%
    filter(Dox == "A") %>%
    select(batch, Cell_line, subnumber, Primer, two_fold_A = two_fold)

  dox_b <- targets %>%
    filter(Dox == "B") %>%
    select(batch, Cell_line, subnumber, Primer, two_fold_B = two_fold)

  ratio_lookup <- full_join(dox_a, dox_b, by = c("batch", "Cell_line", "subnumber", "Primer")) %>%
    mutate(ratio = two_fold_B / two_fold_A) %>%
    select(batch, Cell_line, subnumber, Primer, ratio)

  targets %>%
    left_join(ratio_lookup, by = c("batch", "Cell_line", "subnumber", "Primer")) %>%
    mutate(ratio = ifelse(Dox == "A", ratio, NA))
}

# ---- Main bulk processing function ----
process_qpcr_files <- function(
    file_paths,
    cell_line_filter = NULL,
    output_folder = ".",
    two_fold_title = "2-Fold Relative Increase by Dox Condition, Cell Line and Primer",
    two_fold_subtitle = NULL,
    two_fold_filename = "two_fold_graph.png",
    ratio_title = "Dox+/Dox- Ratio by Subnumber, Cell Line and Primer",
    ratio_subtitle = NULL,
    ratio_filename = "ratio_graph.png",
    ratio_break_point = 10,
    ratio_squish_factor = 6,
    ratio_reference_line = 5,
    ratio_ylim_max = 50,
    box_title = "Dox+/Dox- Ratio Distribution by Cell Line and Primer",
    box_subtitle = NULL,
    box_filename = "ratio_boxplot.png",
    box_break_point = 15,
    box_upper_label = 40,
    box_squish_factor = 6,
    box_reference_line = 5,
    box_ylim_max = 50
) {

  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  two_fold_filename <- file.path(output_folder, two_fold_filename)
  ratio_filename    <- file.path(output_folder, ratio_filename)
  box_filename      <- file.path(output_folder, box_filename)

  raw_data <- bind_rows(lapply(file_paths, read_qpcr_csv))

  if (!is.null(cell_line_filter) && length(cell_line_filter) > 0) {
    raw_data <- raw_data %>% filter(Cell_line %in% as.character(cell_line_filter))
  }

  all_data <- calculate_ddct(raw_data)

  # Explicit subnumber numerical ordering
  all_data <- all_data %>%
    mutate(subnumber_num = suppressWarnings(as.numeric(subnumber))) %>%
    arrange(subnumber_num) %>%
    mutate(subnumber = factor(subnumber, levels = unique(subnumber[order(subnumber_num)])))

  # Group related Primers adjacent to each other
  primer_adjacent_order <- c(
    "17F_16R", "Nfix_17F_16R",
    "17F_18R", "Nfix_17F_18R",
    "17F_sp6",
    "MP_2F_4R", 
    "15F_flag1R", "17F_flag1R", "Nfix_15F_flag1R", "Nfix_17F_flag1R"
  )
  
  present_primers <- unique(all_data$Primer)
  ordered_primers <- c(intersect(primer_adjacent_order, present_primers), setdiff(present_primers, primer_adjacent_order))
  all_data$Primer  <- factor(all_data$Primer, levels = ordered_primers)

  all_data$Dox_condition <- ifelse(all_data$Dox == "A", "-", "+")
  all_data$Cell_line <- as.character(all_data$Cell_line)

  two_fold_upper <- quantile(all_data$two_fold, 0.95, na.rm = TRUE)

  # Graph 1: Two-Fold Plot
  p1 <- ggplot(
    all_data %>% filter(!is.na(two_fold)),
    aes(x = Dox_condition, y = two_fold, fill = Dox_condition)
  ) +
    geom_boxplot() +
    facet_wrap(Primer ~ Cell_line) +
    coord_cartesian(ylim = c(0, two_fold_upper)) +
    theme_minimal() +
    labs(
      title = two_fold_title,
      subtitle = two_fold_subtitle,
      x = "Dox Condition",
      y = "2-Fold Relative Increase (2^-dCt vs Gapdh)"
    ) +
    theme(legend.position = "none")

  ggsave(two_fold_filename, p1, width = 10, height = 7)

  # Consolidated ratios
  ratios <- all_data %>%
    filter(!is.na(ratio)) %>%
    select(Cell_line, subnumber, Primer, ratio) %>%
    distinct()

  ratios$Primer <- factor(ratios$Primer, levels = ordered_primers)

  # Graph 2: Barplot faceted by Primer + Cell Line
  squish_trans <- make_squish_trans(ratio_break_point, ratio_squish_factor)
  upper_max <- suppressWarnings(max(ratios$ratio, na.rm = TRUE))
  label_ceiling <- if (!is.na(ratio_ylim_max)) min(upper_max, ratio_ylim_max, na.rm = TRUE) else upper_max
  upper_breaks <- if (is.finite(label_ceiling) && label_ceiling > ratio_break_point) {
    candidates <- pretty(c(ratio_break_point, label_ceiling), n = 4)
    candidates[candidates > ratio_break_point & candidates <= label_ceiling]
  } else {
    numeric(0)
  }
  y_breaks <- sort(unique(c(0, 5, ratio_break_point, upper_breaks)))

  p2 <- ggplot(ratios, aes(x = subnumber, y = ratio, fill = Primer)) +
    geom_bar(stat = "identity", position = "dodge")

  if (!is.na(ratio_reference_line)) {
    p2 <- p2 + geom_hline(yintercept = ratio_reference_line, linetype = "dotted", linewidth = 0.6, color = "black")
  }

  p2 <- p2 +
    facet_wrap(~ Primer + Cell_line, scales = "free_x") +
    scale_y_continuous(trans = squish_trans, breaks = y_breaks, expand = expansion(mult = c(0.02, 0))) +
    theme_minimal() +
    labs(
      title = ratio_title,
      subtitle = ratio_subtitle,
      x = "Subnumber (Ascending Order)",
      y = "Dox+/Dox- (ddCt fold change)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.na(ratio_ylim_max)) {
    p2 <- p2 + coord_cartesian(ylim = c(0, ratio_ylim_max))
  }

  ggsave(ratio_filename, p2, width = 12, height = 7)

  # ---- Safe Statistical Calculations for Graph 3 ----
  # Calculate 1-Sample t-test vs baseline (mu = 1) for each Primer + Cell_line group
  p_vals_summary <- ratios %>%
    group_by(Primer, Cell_line) %>%
    summarise(
      n = n(),
      p_val = if(n > 1) suppressWarnings(t.test(ratio, mu = 1)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p_val) ~ "",
        p_val < 0.001 ~ "p < 0.001",
        p_val < 0.05  ~ sprintf("p = %.3f", p_val),
        TRUE          ~ sprintf("p = %.2f", p_val)
      )
    )

  # Calculate median values for positioning labels outside boxes
  median_summary <- ratios %>%
    group_by(Primer, Cell_line) %>%
    summarise(
      median_val = median(ratio, na.rm = TRUE),
      .groups = "drop"
    )

  # Graph 3: Box Plot with Nudged Medians & Cell Line Aligned p-Values
  box_trans <- make_squish_trans(box_break_point, box_squish_factor)
  box_breaks <- sort(unique(c(0, 5, 10, box_break_point, box_upper_label)))

  p3 <- ggplot(ratios, aes(x = Cell_line, y = ratio, fill = Cell_line)) +
    # Whisker end-caps (T-bars)
    stat_boxplot(geom = "errorbar", width = 0.2, color = "black", linewidth = 0.6) +
    # Main Boxplot (width = 0.4, spans x - 0.20 to x + 0.20)
    geom_boxplot(outlier.shape = NA, alpha = 0.6, color = "black", width = 0.4) +
    # Individual points jittered
    geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.5, color = "black") +
    # Median labels shifted OUTSIDE the box (nudge_x = 0.28 pushes past the 0.20 right edge)
    geom_text(
      data = median_summary,
      aes(x = Cell_line, y = median_val, label = sprintf("%.2f", median_val)),
      nudge_x = 0.28,
      hjust = 0, # Left-aligned so label starts directly after nudge point
      vjust = 0.5,
      size = 3.8,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    ) +
    # p-values aligned above each respective Cell_line along the x-axis
    geom_text(
      data = p_vals_summary,
      aes(x = Cell_line, y = Inf, label = p_label),
      vjust = 1.5, # Pulls text down slightly from upper margin boundary
      size = 3.6,
      fontface = "italic",
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_wrap(~ Primer, scales = "free_x") +
    scale_y_continuous(trans = box_trans, breaks = box_breaks, expand = expansion(mult = c(0.02, 0.1))) +
    theme_minimal() +
    labs(
      title = box_title,
      subtitle = "Values right of boxes = Medians. Top labels = p-value vs baseline (fold change = 1).",
      x = "Cell Line",
      y = "Dox+/Dox- (ddCt fold change)"
    ) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.na(box_ylim_max)) {
    p3 <- p3 + coord_cartesian(ylim = c(0, box_ylim_max))
  }

  ggsave(box_filename, p3, width = 12, height = 7)

  # Excel Output Generation
  summary_stats <- ratios %>%
    group_by(Cell_line, Primer) %>%
    summarise(
      n = n(),
      mean_ratio = mean(ratio, na.rm = TRUE),
      median_ratio = median(ratio, na.rm = TRUE),
      sd_ratio = sd(ratio, na.rm = TRUE),
      sem_ratio = sd_ratio / sqrt(n),
      .groups = "drop"
    )

  output_file <- file.path(output_folder, "Nfix_MP_bulk_qPCR_output.xlsx")
  wb <- createWorkbook()

  addWorksheet(wb, "Processed_Data")
  writeData(wb, "Processed_Data", all_data)

  addWorksheet(wb, "Ratios")
  writeData(wb, "Ratios", ratios)

  addWorksheet(wb, "Stats_Ratios_Summary")
  writeData(wb, "Stats_Ratios_Summary", summary_stats, startRow = 1)

  addWorksheet(wb, "Graphs")
  insertImage(wb, "Graphs", two_fold_filename, startRow = 1, startCol = 1, width = 10, height = 7, units = "in")
  insertImage(wb, "Graphs", ratio_filename, startRow = 22, startCol = 1, width = 12, height = 7, units = "in")
  insertImage(wb, "Graphs", box_filename, startRow = 43, startCol = 1, width = 12, height = 7, units = "in")

  saveWorkbook(wb, output_file, overwrite = TRUE)

  cat("Processing complete. Median values nudged outside box and p-values aligned per cell line.\n")
}

# Run script
input_folder <- "."
csv_files <- list.files(input_folder, pattern = "_results\\.csv$", full.names = TRUE)

process_qpcr_files(csv_files, output_folder = input_folder)