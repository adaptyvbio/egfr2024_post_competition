#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)

font_add_google("Roboto", "Roboto")
showtext_auto()

# Define colors for secondary structures
dssp_colors <- c(
  "helix" = "#8CD2F4",
  "sheet" = "#56A6D4",
  "loop" = "#8B90DD"
)

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(color = "black"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
}

option_list <- list(
  make_option(c("--input"), type="character", default="/Users/tudorcotet/Documents/Adaptyv/competition_paper/tt.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", default="./plots/dssp",
              help="Output directory path"),
  make_option(c("--type"), type="character", default="both",
              help="Type of DSSP analysis: 'interface', 'whole', or 'both'"),
  make_option(c("--binding"), type="character", default="all",
              help="Filter by binding (Yes, No, or all)"),
  make_option(c("--binding_strength"), type="character", default="all",
              help="Filter by binding strength (Weak, Medium, Strong, or all)"),
  make_option(c("--format"), type="character", default="png",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=3600,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=3600,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI)"),
  make_option(c("--round"), type="character", default="both",
              help="Filter by round (1, 2, or both)"),
  make_option(c("--title"), type="character", default=NULL,
              help="Custom plot title"),
  make_option(c("--subtitle"), type="character", default=NULL,
              help="Custom plot subtitle")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

main <- function() {
  # Read data
  df <- fread(opt$input)
  
  # Filter by round if specified
  if (opt$round != "both") {
    round_num <- as.numeric(opt$round)
    df <- df[round == round_num]
  }
  
  # Filter by binding if specified
  if (opt$binding != "all") {
    df <- df[binding == opt$binding]
  }
  
  # Filter by binding strength if specified
  if (opt$binding_strength != "all") {
    df <- df[binding_strength == opt$binding_strength]
  }
  
  # Prepare DSSP data based on type
  if (opt$type == "both") {
    # Prepare data for both interface and whole structure
    interface_data <- data.table(
      structure = c("helix", "sheet", "loop"),
      percentage = c(
        mean(df$i_helix_pct, na.rm = TRUE),
        mean(df$i_sheet_pct, na.rm = TRUE),
        mean(df$i_loop_pct, na.rm = TRUE)
      ),
      type = "Interface"
    )
    
    whole_data <- data.table(
      structure = c("helix", "sheet", "loop"),
      percentage = c(
        mean(df$helix_pct, na.rm = TRUE),
        mean(df$sheet_pct, na.rm = TRUE),
        mean(df$loop_pct, na.rm = TRUE)
      ),
      type = "Overall"
    )
    
    # Normalize percentages for each type
    interface_data[, percentage := (percentage / sum(percentage)) * 100]
    whole_data[, percentage := (percentage / sum(percentage)) * 100]
    
    # Combine data
    dssp_data <- rbind(interface_data, whole_data)
    
    # Add positions for labels
    dssp_data[type == "Interface", x := 1.5]
    dssp_data[type == "Overall", x := 2.5]
    dssp_data[, ypos := cumsum(percentage) - 0.5 * percentage, by = type]
    
    # Create nested donut plot
    p <- ggplot(dssp_data, aes(x = x, y = percentage, fill = structure)) +
      geom_bar(stat = "identity", width = 0.8, color = "white", linewidth = 0.5) +
      geom_text(aes(y = ypos,
                    label = sprintf("%s: %.1f%%", structure, percentage)),
                color = "black",
                size = 4,
                family = "Roboto",
                position = position_nudge(x = 0.3)) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = dssp_colors,
                       name = "Structure type") +
      scale_x_continuous(limits = c(0, 3.5)) +
      labs(title = "Secondary structure distribution",
           subtitle = sprintf("Round %s - Interface (inner) and Overall (outer)", opt$round)) +
      adaptyv_theme()
    
  } else {
    # Original single donut plot code
    if (opt$type == "interface") {
      dssp_data <- data.table(
        structure = c("helix", "sheet", "loop"),
        percentage = c(
          mean(df$i_helix_pct, na.rm = TRUE),
          mean(df$i_sheet_pct, na.rm = TRUE),
          mean(df$i_loop_pct, na.rm = TRUE)
        )
      )
    } else {
      dssp_data <- data.table(
        structure = c("helix", "sheet", "loop"),
        percentage = c(
          mean(df$helix_pct, na.rm = TRUE),
          mean(df$sheet_pct, na.rm = TRUE),
          mean(df$loop_pct, na.rm = TRUE)
        )
      )
    }
    
    # Normalize percentages
    dssp_data[, percentage := (percentage / sum(percentage)) * 100]
    
    # Add position for labels
    dssp_data[, ypos := cumsum(percentage) - 0.5 * percentage]
    
    # Set title and subtitle
    title <- ifelse(is.null(opt$title),
                   sprintf("%s secondary structure distribution", 
                          ifelse(opt$type == "interface", "Interface", "Overall")),
                   opt$title)
    
    subtitle <- ifelse(is.null(opt$subtitle),
                      sprintf("Round %s", opt$round),
                      opt$subtitle)
    
    # Create single donut plot
    p <- ggplot(dssp_data, aes(x = 2, y = percentage, fill = structure)) +
      geom_bar(stat = "identity", width = 0.8, color = "white", linewidth = 0.5) +
      geom_text(aes(y = ypos,
                    label = sprintf("%s: %.1f%%", structure, percentage)),
                color = "black",
                size = 4,
                family = "Roboto",
                position = position_nudge(x = 0.3)) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = dssp_colors,
                       name = "Structure type") +
      xlim(0.8, 2.5) +
      labs(title = title,
           subtitle = subtitle) +
      adaptyv_theme()
  }
  
  # Create output directory if needed
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  
  # Save plot
  filename <- file.path(opt$output, 
                       sprintf("dssp_%s_round%s_binding%s_strength%s.%s", 
                             opt$type, opt$round, opt$binding, opt$binding_strength, opt$format))
  
  if (opt$format == "png") {
    png(filename, width = opt$width, height = opt$height, res = opt$res)
  } else {
    svg(filename)
  }
  print(p)
  dev.off()
}

if (sys.nframe() == 0) {
  main()
}
