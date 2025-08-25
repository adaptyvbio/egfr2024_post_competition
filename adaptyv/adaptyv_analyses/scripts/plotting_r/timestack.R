#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(ggrepel)  # For better label placement
library(tidyr)    # For replace_na function

# Reuse existing theme and colors from other scripts
font_add_google("Roboto", "Roboto")
showtext_auto()

# Use the exact same theme as barplots.R but remove grid lines
adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 16, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5),  # Horizontal x-axis labels
      axis.text.y = element_text(size = 12),                          # Increased y-axis text size
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),                                   # Remove all grid lines
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Define color palettes (reused from existing scripts)
binding_colors <- c(
  "None" = "#E5E7EB",
  "Not expressed" = "#F9F5C9", 
  "Weak" = "#C4B2FB",
  "Medium" = "#A7C1FB",
  "Strong" = "#2A4DD0",
  "Missing binding data" = "#3E6175"
)

expression_colors <- c(
  "None" = '#E5E7EB',
  "Low" = "#9EA2AF",
  "Medium" = "#E9C435", 
  "High" = "#69B7A7"
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

# Add command line options for filtering
option_list <- list(
  make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", default="./plots/timestacks/",
              help="Output directory path"),
  make_option(c("--category_column"), type="character", default="design_category",
              help="Column to use for categories (binding_strength, expression, selected)"),
  make_option(c("--round"), type="character", default="both",
              help="Filter by round (1, 2, or both)"),
  make_option(c("--binding"), type="character", default="all",
              help="Filter by binding status (Yes, No, or all)"),
  make_option(c("--binding_strength"), type="character", default="all", 
              help="Filter by binding strength (Strong, Medium, Weak, or all)"),
  make_option(c("--expression"), type="character", default="High",
              help="Filter by expression level (High, Medium, Low, or all)"),
  make_option(c("--de_novo"), type="character", default="all",
              help="Filter by de novo status (De novo, Existing binder, or all)"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=1800,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=1400,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=300,
              help="Plot resolution (DPI)"),
  make_option(c("--title"), type="character", default="Highly expressed designs",
              help="Custom plot title"),
  make_option(c("--subtitle"), type="character", default="Evolution of design categories across rounds",
              help="Custom plot subtitle"),
  make_option(c("--remove_not_mentioned"), type="logical", default=TRUE,
              help="Remove 'Not mentioned' category"),
  make_option(c("--max_categories"), type="integer", default=0,
              help="Maximum number of categories to show (0 for all)"),
  make_option(c("--reference_round"), type="integer", default=1,
              help="Round to use for selecting top categories")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

main <- function() {
  tryCatch({
    log_info('Starting script execution')
    
    # Read data
    log_info(sprintf('Loading data from %s', opt$input))
    data <- fread(opt$input)
    
    # Validate category column exists
    if (!(opt$category_column %in% names(data))) {
      stop(sprintf("Category column '%s' not found in data", opt$category_column))
    }
    
    # Filter by round if specified
    if (opt$round != "both") {
      round_num <- as.numeric(opt$round)
      data <- data[round == round_num]
      log_info(sprintf('Filtered data for round %d', round_num))
    }
    
    # Apply filters
    if (opt$binding != "all") {
      data <- data[binding == opt$binding]
      log_info(sprintf('Filtered for binding status: %s', opt$binding))
    }
    
    if (opt$binding_strength != "all") {
      data <- data[binding_strength == opt$binding_strength]
      log_info(sprintf('Filtered for binding strength: %s', opt$binding_strength))
    }
    
    if (opt$expression != "all") {
      data <- data[expression == opt$expression]
      log_info(sprintf('Filtered for expression level: %s', opt$expression))
    }
    
    if (opt$de_novo != "all") {
      data <- data[de_novo == opt$de_novo]
      log_info(sprintf('Filtered for de novo status: %s', opt$de_novo))
    }
    
    # Remove "Not mentioned" if specified
    if (opt$remove_not_mentioned) {
      data <- data[get(opt$category_column) != "Not mentioned"]
      log_info('Removed "Not mentioned" category')
    }
    
    # Calculate percentages per round using dplyr
    plot_data <- data %>%
      group_by(round, !!sym(opt$category_column)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(round) %>%
      mutate(percentage = count / sum(count) * 100) %>%
      ungroup()
    
    # Get all unique categories and rounds
    all_categories <- unique(plot_data[[opt$category_column]])
    all_rounds <- unique(plot_data$round)
    
    # Create complete grid with 0% for missing combinations
    plot_data <- expand.grid(
      round = all_rounds,
      category = all_categories
    ) %>%
      setNames(c("round", opt$category_column)) %>%
      left_join(plot_data, by = c("round", opt$category_column)) %>%
      mutate(
        count = ifelse(is.na(count), 0, count),
        percentage = ifelse(is.na(percentage), 0, percentage)
      )
    
    # Convert round to factor
    plot_data$round <- factor(plot_data$round)
    
    # Select appropriate color palette based on category
    color_values <- switch(opt$category_column,
                          "binding_strength" = binding_colors,
                          "expression" = expression_colors,
                          "selected" = selection_status_colors,
                          "design_category" = design_category_colors,
                          setNames(rainbow(length(unique(plot_data[[opt$category_column]]))), 
                                 unique(plot_data[[opt$category_column]])))  # Fallback color palette
    
    # Ensure color_values only contains colors for categories that exist in the filtered data
    color_values <- color_values[names(color_values) %in% unique(plot_data[[opt$category_column]])]
    
    # Set factor levels based on category
    category_levels <- switch(opt$category_column,
      "binding_strength" = c("None", "Not expressed", "Weak", "Medium", "Strong"),
      "expression" = c("None", "Low", "Medium", "High"),
      "selected" = c("Top 100", "Adaptyv selection", "Not selected"),
      "design_category" = c("De novo", "Optimized binder", "Diversified binder", "Hallucination", "Physics-based")
    )
    
    # Filter categories if max_categories is set
    if (opt$max_categories > 0) {
      # Get top categories based on percentage in reference round
      top_categories <- plot_data %>%
        filter(round == opt$reference_round) %>%
        arrange(desc(percentage)) %>%
        head(opt$max_categories) %>%
        pull(!!sym(opt$category_column))
      
      # Filter plot_data to keep only top categories
      plot_data <- plot_data %>%
        filter(!!sym(opt$category_column) %in% top_categories)
      
      log_info(sprintf('Filtered to show top %d categories based on round %d', 
                      opt$max_categories, opt$reference_round))
    }
    
    # Set factor levels based on category (modified section)
    if (!is.null(category_levels)) {
      # Get the actual levels present in the data
      existing_levels <- unique(plot_data[[opt$category_column]])
      # Only use the levels that exist in the data
      valid_levels <- category_levels[category_levels %in% existing_levels]
      # If max_categories is set, order by percentage in reference round
      if (opt$max_categories > 0) {
        level_order <- plot_data %>%
          filter(round == opt$reference_round) %>%
          arrange(desc(percentage)) %>%
          pull(!!sym(opt$category_column))
        valid_levels <- valid_levels[valid_levels %in% level_order]
      }
      # Set the factor levels
      plot_data[[opt$category_column]] <- factor(plot_data[[opt$category_column]], 
                                               levels = valid_levels)
    }
    
    # Set default title if not provided
    if (opt$title == "") {
      title <- sprintf("%s distribution", stringr::str_to_sentence(gsub("_", " ", opt$category_column)))
    } else {
      title <- opt$title
    }

    if (opt$subtitle == "") {
      subtitle <- sprintf("Evolution across rounds")
    } else {
      subtitle <- opt$subtitle
    }
    
    # Create slope chart
    p <- ggplot(plot_data, aes(x = round, 
                              y = percentage, 
                              group = !!sym(opt$category_column),
                              color = !!sym(opt$category_column))) +
      # Add connecting lines
      geom_line(linewidth = 1.5, alpha = 0.8) +
      # Add points
      geom_point(size = 4, alpha = 0.9) +
      # Add percentage labels
      geom_text_repel(
        aes(label = sprintf("%.1f%%", percentage)),
        size = 4,
        fontface = "bold",
        box.padding = 0.8,
        point.padding = 0.5,
        min.segment.length = 0.2,
        segment.color = "grey50",
        segment.size = 0.3,
        max.overlaps = Inf,
        show.legend = FALSE,
        direction = "y",
        nudge_x = 0.1
      ) +
      # Enhanced color styling
      scale_color_manual(values = color_values,
                        name = stringr::str_to_sentence(gsub("_", " ", opt$category_column))) +
      # Improved labels
      labs(title = title,
           subtitle = subtitle,
           x = "Round",
           y = "Percentage (%)") +
      # Use the consistent theme
      adaptyv_theme() +
      # Scale adjustments
      scale_y_continuous(
        limits = function(x) c(0, max(x) * 1.1),
        expand = expansion(mult = c(0, 0.1))
      ) +
      scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
      # Override specific theme elements for this plot
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),  # Keep x-axis labels horizontal
        panel.grid.major.x = element_blank()                  # Remove vertical gridlines
      )

    # Create output directory if it doesn't exist
    dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
    
    # Generate filename
    filename <- file.path(opt$output, 
                         sprintf("timestack_%s.%s", 
                                opt$category_column, 
                                opt$format))
    
    # Save plot
    if (opt$format == "png") {
      png(filename, width = opt$width, height = opt$height, res = opt$res)
      print(p)
      dev.off()
      log_info(sprintf('Saved PNG plot to %s', filename))
    } else if (opt$format == "svg") {
      svg(filename, width = opt$width / opt$res, height = opt$height / opt$res)
      print(p)
      dev.off()
      log_info(sprintf('Saved SVG plot to %s', filename))
    }
    
    log_info('Script completed successfully')
    
  }, error = function(e) {
    log_error(sprintf("Script failed with error: %s", e$message))
    stop(e)
  })
}

if (sys.nframe() == 0) {
  main()
}