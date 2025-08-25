#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)

# Import theme and colors
font_add_google("Roboto", "Roboto")
showtext_auto()

# Format metric names consistently with other plots
format_metric_name <- function(metric, for_title = FALSE) {
    metric_formats <- list(
        "kd" = if(for_title) "KD" else "'K'[D]",
        "-log10_kd" = if(for_title) "-log10(KD)" else "-log[10]('K'[D])",
        "normalized_kd" = if(for_title) "Normalized KD" else "paste('Normalized ', 'K'[D])",
        "iptm" = "ipTM",
        "pae_interaction" = "iPAE",
        "esm_pll" = "ESM PLL",
        "normalized_esm_pll" = "Normalized ESM PLL",
        "plddt" = "pLDDT",
        "pae_score" = "PAE Score",
        "sequence_length" = "Sequence length",
        "interface_nres" = "Interface residues"
    )
    
    metric_lower <- tolower(metric)
    if (metric_lower %in% names(metric_formats)) {
        formatted <- metric_formats[[metric_lower]]
        return(if(!for_title && !grepl("^-log", metric_lower)) sprintf("'%s'", formatted) else formatted)
    } else {
        return(stringr::str_to_title(gsub("_", " ", metric)))
    }
}

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Command line options
option_list <- list(
  make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", default="./plots/densities/",
              help="Output directory path"),
  make_option(c("--metric"), type="character", default="sequence_length",
              help="Metric to plot (e.g., iptm, pae_interaction, esm_pll)"),
  make_option(c("--category"), type="character", default="round",
              help="Category to split by (e.g., binding, design_category, round)"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=1600,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=1200,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=300,
              help="Plot resolution (DPI)"),    
  make_option(c("--round"), type="character", default='both',
              help="Filter by round (1, 2, or both)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Define color palettes
binding_colors <- c(
  "None" = "#E5E7EB",
  "Unknown" = "#F9F5C9", 
  "Weak" = "#C4B2FB",
  "Medium" = "#A7C1FB",
  "Strong" = "#2A4DD0",
  "Missing binding data" = "#3E6175",
  "% binders" = "#DC7A73"
)

expression_colors <- c(
  "Low" = "#9EA2AF",
  "Medium" = "#E9C435", 
  "High" = "#69B7A7",
  "Not expressed" = '#DC7A73'
)

selection_status_colors <- c(
  "Top 100" = "#8B90DD",
  "Adaptyv selection" = "#8CD2F4",
  "Not selected" = "#3E6175"
)

design_category_colors <- c(
  "De novo" = "#56A6D4",
  "Optimized binder" = "#8B90DD",
  "Diversified binder" = "#8CD2F4",
  "Hallucination" = "#3E6175",
  "Physics-based" = "#69B7A7"
)

# Add round colors
round_colors <- c(
    "2" = "#8B90DD",  # purple
    "1" = "#8CD2F4"   # blue
)

binding_status_colors <- c(
  "Binder" = "#56A6D4",
  "Non-binder" = "#E5E7EB"
)

adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, 
                                 selection_status_colors, design_category_colors)))

main <- function() {
  tryCatch({
    # Read data
    df <- fread(opt$input)
    
    # Filter by round if specified
    if (opt$round != "both") {
      round_num <- as.numeric(opt$round)
      df <- df[round == round_num]
      log_info(sprintf('Filtered data for round %d', round_num))
    }
    
    # Create category based on input parameter
    if (opt$category == "binding") {
      df$plot_category <- ifelse(df$binding == "Yes", "Binder", "Non-binder")
      color_values <- binding_status_colors
      category_title <- "Binding status"
    } else if (opt$category == "design_category") {
      df$plot_category <- df$design_category
      color_values <- design_category_colors
      category_title <- "Design category"
    } else if (opt$category == "round") {
      df$plot_category <- as.character(df$round)
      color_values <- round_colors
      category_title <- "Round"
    } else {
      df$plot_category <- df[[opt$category]]
      color_values <- adaptyv_colors
      category_title <- stringr::str_to_title(gsub("_", " ", opt$category))
    }
    
    # Remove NA categories
    df <- df[!is.na(plot_category)]
    
    # Format metric name for title
    metric_title <- format_metric_name(opt$metric, for_title = TRUE)
    
    # Create overlapping density plot
    p <- ggplot(df, aes(x = get(opt$metric), fill = plot_category)) +
      geom_density(alpha = 0.4) +
      geom_density(aes(color = plot_category), fill = NA, linewidth = 0.8) +
      scale_fill_manual(values = color_values) +
      scale_color_manual(values = color_values) +
      labs(
        title = sprintf("%s", metric_title),
        subtitle = sprintf("Shift across rounds"),
        x = metric_title,
        y = "Density",
        fill = category_title,
        color = category_title
      ) +
      adaptyv_theme() +
      guides(color = "none")
    
    # Create output directory if needed
    dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
    
    # Save plot
    filename <- file.path(opt$output, 
                         sprintf("density_%s_by_%s_round%s.%s", 
                                opt$metric,
                                opt$category,
                                opt$round,
                                opt$format))
    
    if (opt$format == "png") {
      png(filename, width = opt$width, height = opt$height, res = opt$res)
    } else {
      svg(filename, width = opt$width / opt$res, height = opt$height / opt$res)
    }
    print(p)
    dev.off()
    
    log_info(sprintf("Plot saved to %s", filename))
    
  }, error = function(e) {
    log_error(sprintf("Error: %s", e$message))
    stop(e)
  })
}

if (sys.nframe() == 0) {
  main()
}