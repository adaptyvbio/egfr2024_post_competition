#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)

# Add Roboto font
font_add_google("Roboto", "Roboto")
showtext_auto()

# Command line options
option_list <- list(
  make_option(c("--input"), type="character", 
              default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", 
              default="./plots/binding_affinity/",
              help="Output directory path"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=6000,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=4000,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Define colors for rounds
round_colors <- c(
  "1" = "#67AFD8",  # Light blue
  "2" = "#8D92DE"   # Purple
)

# Theme for the plot
adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 16, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.position = c(1.05, 1.2),  # Moved legend to x=1.05, y=1.2
      legend.justification = c(1, 1),    # Top right anchor
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(30, 50, 60, 30)  # Increased right margin to accommodate legend
    )
}

main <- function() {
  tryCatch({
    # Read data
    log_info(sprintf('Loading data from %s', opt$input))
    data <- fread(opt$input)
    
    # Filter for sequences with KD values and sort
    plot_data <- data[!is.na(kd)]
    plot_data <- plot_data[order(kd)][1:60]

    # Get the range of KD values
    kd_range <- range(plot_data$kd)
    y_min <- floor(log10(kd_range[1]))  # Round down to nearest power of 10
    y_max <- ceiling(log10(kd_range[2]))  # Round up to nearest power of 10
    
    # Add sequence index for x-axis (not reversed anymore)
    plot_data[, sequence_index := seq_len(.N)]
    
    # Trim sequence names to 20 characters
    plot_data[, name := ifelse(nchar(name) > 20, 
                              paste0(substr(name, 1, 20), "..."), 
                              name)]
    
    # Count total sequences and characterized ones
    total_sequences <- nrow(data)
    characterized_sequences <- nrow(data[!is.na(kd)])
    
    # Control KD values from violin_plots.R
    cetuximab_kd <- mean(c(9.94e-9, 3.33e-9))  # ~6.635e-9
    egf_kd <- mean(c(7.95e-7, 7.57e-7, 7.25e-7))  # ~7.59e-7
    
    # Reference KD values for common antibody-antigen interactions
    reference_kds <- data.table(
      value = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5),
      label = c("(1 nM)", "(10 nM)", "(100 nM)", "(1 µM)", "(10 µM)")
    )
    
    # Create plot
    p <- ggplot(plot_data) +
      # Add background rectangles for binding regions
      annotate("rect", xmin = -Inf, xmax = Inf, 
               ymin = 1e-5, ymax = 3e-5,  # Adjusted to match data
               fill = "#F8F0FF", alpha = 0.3) +
      # Add weak binding threshold line at the bottom
      geom_hline(yintercept = 1e-5,
                 color = "#E9ECEF",
                 linetype = "dashed") +
      # Add reference KD lines and labels
      geom_hline(data = reference_kds,
                 aes(yintercept = value),
                 color = "#E9ECEF",
                 linetype = "dotted") +
      geom_text(data = reference_kds,
                aes(x = 0, y = value, label = label),
                hjust = 0, vjust = -0.5,
                color = "#6B7280", size = 2.5) +
      # Add control KD lines and labels
      geom_hline(yintercept = cetuximab_kd,
                 linetype = "dashed",
                 color = "#E9C435",
                 linewidth = 0.5) +
      geom_hline(yintercept = egf_kd,
                 linetype = "dashed",
                 color = "#69B7A7",
                 linewidth = 0.5) +
      annotate("text",
               x = Inf,
               y = cetuximab_kd,
               label = "Cetuximab control",
               hjust = 1.1,
               vjust = -0.5,
               size = 3,
               color = "#E9C435") +
      annotate("text",
               x = Inf,
               y = egf_kd,
               label = "EGF control",
               hjust = 1.1,
               vjust = -0.5,
               size = 3,
               color = "#69B7A7") +
      # Add binding threshold labels
      annotate("text", x = 2, y = 2e-5, label = "Weak binding", 
               hjust = 0, vjust = 0.5, color = "#A4A7B7", size = 3) +
      # Customize scales
      scale_x_continuous(
        breaks = plot_data$sequence_index,
        labels = plot_data$name,
        expand = expansion(mult = c(0.02, 0.02))  # Minimal padding
      ) +
      scale_y_log10(
        breaks = 10^seq(y_min, y_max),
        labels = parse(text = sprintf("10^%d", seq(y_min, y_max))),
        expand = expansion(mult = c(0.02, 0.02))  # Minimal padding
      ) +
      # Use coord_cartesian for limits to prevent data clipping
      coord_cartesian(
        ylim = c(3e-5, 1e-9),  # Reversed limits for better visualization
        clip = "off"  # Allow annotations outside plot area
      ) +
      scale_color_manual(values = round_colors, name = "Round") +
      # Add labels
      labs(
        title = "Binding affinity",
        subtitle = sprintf("All binders out of 601 / %d sequences characterized", 
                          total_sequences),
        x = NULL,
        y = "KD [M]"
      ) +
      # Add theme
      adaptyv_theme() +
      # Add points for each measurement (moved to end to be on top)
      geom_point(aes(x = sequence_index, y = kd, color = factor(round)),
                size = 1.5)
    
    # Create output directory
    dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
    
    # Save plot
    filename <- file.path(opt$output, 
                         sprintf("binding_affinity.%s", opt$format))
    
    ggsave(filename, 
           p, 
           width = opt$width/opt$res, 
           height = opt$height/opt$res, 
           dpi = opt$res)
    
    log_info(sprintf('Plot saved to %s', filename))
    
  }, error = function(e) {
    log_error(sprintf('Error: %s', e$message))
    stop(e)
  })
}

if (!interactive()) {
  main()
} 