#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)

# Define theme and colors directly since source file is not accessible
font_add_google("Roboto", "Roboto")
showtext_auto()

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 16, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )
}

binding_colors <- c(
  "None" = "#E5E7EB",
  "Not expressed" = "#9EA2AF", 
  "Weak" = "#C4B2FB",
  "Medium" = "#A7C1FB",
  "Strong" = "#2A4DD0",
  "Missing binding data" = "#3E6175",
  "% binders" = "#DC7A73"
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
  "Not selected" = "#E5E7EB"
)

design_category_colors <- c(
  "De novo" = "#56A6D4",
  "Optimized binder" = "#8B90DD",
  "Diversified binder" = "#8CD2F4",
  "Hallucination" = "#3E6175"
)

# Add new colors for the model type columns
model_type_colors <- c(
  "Yes" = "#56A6D4",
  "No" = "#E5E7EB"
)

# Create specific colors for each model type
bindcraft_colors <- c(
  "BindCraft" = "#56A6D4",
  "Other" = "#E5E7EB"
)

rfdiffusion_colors <- c(
  "RFdiffusion" = "#8B90DD",
  "Other" = "#E5E7EB"
)

esm_colors <- c(
  "ESM" = "#8CD2F4",
  "Other" = "#E5E7EB"
)

# Updated TIMED colors with four categories using design_category_colors
timed_colors <- c(
  "TIMED" = "#56A6D4",         # De novo blue
  "ProteinMPNN" = "#8B90DD",   # Optimized binder purple
  "ESM-IF" = "#8CD2F4",        # Diversified binder light blue
  "Other" = "#E5E7EB"          # Light gray
)

af2_backprop_colors <- c(
  "AF2 backprop" = "#8CD2F4",
  "Other" = "#E5E7EB"
)

other_hallucination_colors <- c(
  "Other hallucination" = "#8CD2F4",
  "Other" = "#E5E7EB"
)

metric_colors <- c("ESM2 PLL" = "#8CD2F4", "ipTM" = "#8B90DD", "iPAE" = "#3E6175")

adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, selection_status_colors, design_category_colors, metric_colors, 
                                 bindcraft_colors, rfdiffusion_colors, esm_colors, timed_colors, af2_backprop_colors, other_hallucination_colors)))

option_list <- list(
  make_option(c("--x_column"), type="character", default='design_category',
              help="Column name to plot on x-axis"),
  make_option(c("--color_column"), type="character", default='selected',
              help="Column name to use for stacked bars coloring"),
  make_option(c("--output"), type="character", default="./plots/barplots/",
              help="Output file path"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg) [default=%default]"),
  make_option(c("--width"), type="integer", default=2600,
              help="Plot width in pixels [default=%default]"),
  make_option(c("--height"), type="integer", default=2200,
              help="Plot height in pixels [default=%default]"),
  make_option(c("--res"), type="integer", default=300,
              help="Plot resolution (DPI) [default=%default]"),
  make_option(c("--input"), type="character", default="./data/processed/all_submissions_new.csv",
              help="Input CSV file path [default=%default]"),
  make_option(c("--remove_not_mentioned"), type="logical", default=FALSE,
              help="Remove 'Not mentioned' category [default=%default]"),
  make_option(c("--remove_missing_data"), type="logical", default=TRUE,
              help="Remove 'Missing data' category [default=%default]"),
  make_option(c("--round"), type="character", default='both',
              help="Filter by round (1, 2, or both) [default=%default]"),
  make_option(c("--title"), type="character", default="Design category",
              help="Custom title for the plot [default=%default]"),
  make_option(c("--top_n"), type="integer", default=10,
              help="Only show top N most abundant categories (0 for all) [default=%default]"),
  make_option(c("--subtitle"), type="character", default="Number of designs submitted vs. selected for validation",
              help="Custom subtitle for the plot [default=%default]"),
  make_option(c("--sort"), type="character", default="size",
              help="Sort x-axis bars by 'size' or 'alpha' [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$x_column) || is.null(opt$color_column) || is.null(opt$output) || is.null(opt$input)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

log_threshold(INFO)
log_formatter(formatter_paste)

main <- function() {
  tryCatch({
    log_info('Starting script execution')
    
    # Check if input file exists
    if (!file.exists(opt$input)) {
      stop(sprintf("Input file does not exist: %s", opt$input))
    }
    
    # Load data
    log_info(sprintf('Loading data file: %s', opt$input))
    df <- read.csv(opt$input)
    
    # Filter by round if specified
    if (opt$round != "both") {
      round_num <- as.numeric(opt$round)
      df <- df[df$round == round_num,]
      log_info(sprintf('Filtered data for round %d', round_num))
    }
    
    # Create plot
    log_info(sprintf('Creating plot with x=%s colored by %s', opt$x_column, opt$color_column))
    
    # Calculate counts and percentages
    plot_data <- as.data.table(df)
    plot_data[is.na(get(opt$color_column)), (opt$color_column) := "missing data"]
    
    # Filter to top N categories if specified
    if (opt$top_n > 0) {
      top_categories <- plot_data[, .N, by = get(opt$x_column)][order(-N)][1:min(opt$top_n, .N)]
      plot_data <- plot_data[get(opt$x_column) %in% top_categories$get]
      log_info(sprintf('Filtered to top %d categories', opt$top_n))
    }
    
    # Remove "Not mentioned" if specified
    if (opt$remove_not_mentioned) {
      plot_data <- plot_data[get(opt$color_column) != "Not mentioned"]
      plot_data <- plot_data[get(opt$x_column) != "Not mentioned"]
      log_info('Removed "Not mentioned" category')
    }

    # Remove "Missing data" if specified
    if (opt$remove_missing_data) {
      plot_data <- plot_data[!is.na(get(opt$color_column)) & get(opt$color_column) != "" & get(opt$color_column) != "Missing data"]
      plot_data <- plot_data[!is.na(get(opt$x_column)) & get(opt$x_column) != "" & get(opt$x_column) != "Missing data"]
      log_info('Removed missing data, empty strings, and NA values')
    }
    
    # Remove "No" from selected column if that's what we're plotting
    if (opt$x_column == "selected") {
      plot_data <- plot_data[get(opt$x_column) != "No"]
    }
    
    # Remove "Unknown" when x axis is "binding"
    if (opt$x_column == "binding") {
      plot_data <- plot_data[get(opt$x_column) != "Unknown"]
      log_info('Removed "Unknown" from binding column')
    }
    
    print(nrow(plot_data))
    counts <- plot_data[, .(
      total = .N,
      binders = sum(get(opt$color_column) %in% c("Weak", "Medium", "Strong")),
      expressed = sum(get(opt$color_column) %in% c("Low", "Medium", "High")),
      selected = if(opt$color_column == "selected") {
        round(100 * sum(get(opt$color_column) %in% c("Top 100", "Adaptyv selection")) / .N)
      } else if(opt$color_column == "binding_strength") {
        round(100 * sum(get(opt$color_column) %in% c("Weak", "Medium", "Strong")) / .N)
      } else if(opt$color_column == "expression") {
        round(100 * sum(get(opt$color_column) %in% c("Low", "Medium", "High")) / .N)
      } else if(opt$color_column == "TIMED") {
        # For the TIMED column, we don't need to calculate a combined percentage
        # Each design is already categorized as either TIMED, ProteinMPNN, ESM-IF, or Other
        # Just to display something consistent, show the percentage of non-"Other" designs
        round(100 * sum(get(opt$color_column) != "Other") / .N)
      } else if(opt$color_column %in% c("BindCraft", "RFdiffusion", "ESM", "AF2_backprop", "Other_hallucination")) {
        # Continue using "Yes" for percentage calculation, since the transformation happens later
        round(100 * sum(get(opt$color_column) == "Yes") / .N)
      } else {
        0
      }
    ), by = get(opt$x_column)]

    setnames(counts, "get", opt$x_column)
    
    # Order categories by total count or alphabetically based on sort argument
    category_order <- if(opt$sort == "size") {
      counts[order(-total)][[opt$x_column]]
    } else if(opt$x_column == "binding") {
      # Order x axis as "Yes" then "No" for binding
      c("Yes", "No")
    } else {
      sort(unique(plot_data[[opt$x_column]]))
    }
    
    plot_data[[opt$x_column]] <- factor(plot_data[[opt$x_column]], levels = category_order)
    counts[[opt$x_column]] <- factor(counts[[opt$x_column]], levels = category_order)
    
    # Select appropriate color palette based on column name
    color_values <- switch(opt$color_column,
                         "binding_strength" = binding_colors,
                         "expression" = expression_colors,
                         "selected" = selection_status_colors,
                         "design_category" = design_category_colors,
                         "BindCraft" = bindcraft_colors,
                         "RFdiffusion" = rfdiffusion_colors,
                         "ESM" = esm_colors,
                         "TIMED" = timed_colors,
                         "AF2_backprop" = af2_backprop_colors,
                         "Other_hallucination" = other_hallucination_colors,
                         adaptyv_colors)
    
    # For binding strength, set specific factor levels
    if (opt$color_column == "binding_strength") {
      plot_data[plot_data[[opt$color_column]] == 'Unknown', opt$color_column] <- 'Not expressed'
      plot_data[[opt$color_column]] <- factor(plot_data[[opt$color_column]], 
                                            levels = c("None", "Not expressed", "Weak", "Medium", "Strong"))
    }
    
    # For expression, set specific factor levels
    if (opt$color_column == "expression") {
      plot_data[[opt$color_column]] <- factor(plot_data[[opt$color_column]], 
                                            levels = c("None", "Low", "Medium", "High"))
    }
    
    # For selected, set specific factor levels
    if (opt$color_column == "selected") {
      plot_data[plot_data[[opt$color_column]] == "No", opt$color_column] <- "Not selected"
      plot_data[[opt$color_column]] <- factor(plot_data[[opt$color_column]], 
                                            levels = c("Top 100", "Adaptyv selection", "Not selected"))
    }
    
    # For binary model type columns, set specific factor levels with descriptive names
    if(opt$color_column == "BindCraft") {
      plot_data$BindCraft[plot_data$BindCraft == "Yes"] <- "BindCraft"
      plot_data$BindCraft[plot_data$BindCraft == "No"] <- "Other"
      plot_data$BindCraft <- factor(plot_data$BindCraft, levels = c("BindCraft", "Other"))
    } else if(opt$color_column == "RFdiffusion") {
      plot_data$RFdiffusion[plot_data$RFdiffusion == "Yes"] <- "RFdiffusion"
      plot_data$RFdiffusion[plot_data$RFdiffusion == "No"] <- "Other"
      plot_data$RFdiffusion <- factor(plot_data$RFdiffusion, levels = c("RFdiffusion", "Other"))
    } else if(opt$color_column == "ESM") {
      plot_data$ESM[plot_data$ESM == "Yes"] <- "ESM"
      plot_data$ESM[plot_data$ESM == "No"] <- "Other"
      plot_data$ESM <- factor(plot_data$ESM, levels = c("ESM", "Other"))
    } else if(opt$color_column == "TIMED") {
      # Handle TIMED with four categories
      # The values should be "TIMED", "ProteinMPNN", "ESM-IF", or "Other" from preprocessing
      # Make sure all values are valid and factorize in the specific order
      valid_categories <- c("TIMED", "ProteinMPNN", "ESM-IF", "Other")
      plot_data$TIMED[!plot_data$TIMED %in% valid_categories] <- "Other"
      plot_data$TIMED <- factor(plot_data$TIMED, levels = valid_categories)
    } else if(opt$color_column == "AF2_backprop") {
      plot_data$AF2_backprop[plot_data$AF2_backprop == "Yes"] <- "AF2 backprop"
      plot_data$AF2_backprop[plot_data$AF2_backprop == "No"] <- "Other"
      plot_data$AF2_backprop <- factor(plot_data$AF2_backprop, levels = c("AF2 backprop", "Other"))
    } else if(opt$color_column == "Other_hallucination") {
      plot_data$Other_hallucination[plot_data$Other_hallucination == "Yes"] <- "Other hallucination"
      plot_data$Other_hallucination[plot_data$Other_hallucination == "No"] <- "Other"
      plot_data$Other_hallucination <- factor(plot_data$Other_hallucination, levels = c("Other hallucination", "Other"))
    }
    
    if (opt$title == "") {
      title <- sprintf("%s representation", stringr::str_to_sentence(gsub("_", " ", opt$x_column)))
    } else {
      title <- opt$title
    }

    if (opt$subtitle == "") {
      subtitle <- sprintf("Colored by %s", gsub("_", " ", opt$color_column))
    } else {
      subtitle <- opt$subtitle
    }

    # Create plot
    p <- ggplot(plot_data) +
      geom_bar(aes(x = get(opt$x_column), fill = get(opt$color_column)), 
              position = "stack", color = "black") +
      geom_text(stat = "count", 
               aes(x = get(opt$x_column), group = get(opt$color_column), 
                   label = after_stat(count)), 
               position = position_stack(vjust = 0.5),
               color = "white",
               size = 5,  # Increased text size from 3 to 5
               check_overlap = TRUE) + # Added check_overlap
      geom_text(data = counts,
               aes(x = get(opt$x_column), y = total, 
                   label = sprintf("%d%% %s", selected, 
                                 ifelse(opt$color_column == "binding_strength", "binders",
                                      ifelse(opt$color_column == "expression", "expressed",
                                           ifelse(opt$color_column == "selected", "selected", 
                                                ifelse(opt$color_column == "TIMED", "design methods",
                                                     ifelse(opt$color_column %in% c("BindCraft", "RFdiffusion", "ESM", "AF2_backprop", "Other_hallucination"),
                                                          paste0(opt$color_column), ""))))))),
               vjust = -0.5,
               color = "#3D7A9A",
               size = 4,
               fontface = "bold") +
      scale_fill_manual(values = color_values,
                      name = ifelse(opt$color_column %in% c("BindCraft", "RFdiffusion", "ESM", "TIMED", "AF2_backprop", "Other_hallucination"),
                                  "Design model",
                                  stringr::str_to_sentence(gsub("_", " ", opt$color_column)))) +
      labs(title = title,
           subtitle = subtitle,
           x = stringr::str_to_sentence(gsub("_", " ", opt$x_column)),
           y = "Number of designs") +
      adaptyv_theme() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      theme(axis.text.x = element_text(size = 14))
    
    # Create output directory if it doesn't exist
    dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
    
    # Generate unique filename based on x and color columns
    filename <- file.path(opt$output, sprintf("barplot_%s_by_%s.%s", 
                                            opt$x_column, opt$color_column, opt$format))
    
    # Save plot
    if (opt$format == "png") {
      png(filename, bg = 'white', 
          height = opt$height, width = opt$width, res = opt$res)
      print(p)
      dev.off()
      log_info(sprintf('Saved PNG plot to %s', filename))
    } else if (opt$format == "svg") {
      svg(filename, bg = 'transparent', width = opt$width/opt$res, height = opt$height/opt$res)
      print(p)
      dev.off()
      log_info(sprintf('Saved SVG plot to %s', filename))
    } else {
      stop(sprintf("Unsupported format: %s", opt$format))
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
