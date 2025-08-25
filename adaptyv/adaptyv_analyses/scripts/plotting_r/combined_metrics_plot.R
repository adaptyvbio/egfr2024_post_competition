#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger, warn.conflicts = FALSE)
library(optparse)
library(gridExtra)
library(grid)

# Import theme and colors
font_add_google("Roboto", "Roboto")
showtext_auto()

# Define color palettes from adaptyv_theme.R
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
  "Not mentioned" = "#F9F5C9"
)

metric_colors <- c("ESM2 PLL" = "#8CD2F4", "ipTM" = "#8B90DD", "iPAE" = "#3E6175")
adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, 
                                 selection_status_colors, design_category_colors, 
                                 metric_colors)))

# Modified theme without gray grid
adaptyv_theme_no_grid <- function() {
    theme_minimal() +
        theme(
            text = element_text(family = "Roboto", color = "#333333"),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
            plot.subtitle = element_markdown(size = 10, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 10)),
            axis.title = element_text(size = 10, face = "bold"),
            axis.title.y = element_text(angle = 90, vjust = 0.5),
            axis.text = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.key = element_rect(fill = "white", color = NA),
            legend.title = element_text(face = "bold", size = 8),
            legend.text = element_text(size = 7),
            legend.position = "bottom",
            legend.box.background = element_rect(color = "white"),
            plot.margin = margin(5, 5, 5, 5)
        )
}

# Define command line options
option_list <- list(
    make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
                help="Input CSV file path"),
    make_option(c("--output"), type="character", default="./plots/combined_metrics/",
                help="Output directory path"),
    make_option(c("--y_column"), type="character", default="kd",
                help="Column name for y-axis (typically KD or another outcome measure)"),
    make_option(c("--color_by"), type="character", default="design_category",
                help="Column name for point colors"),
    make_option(c("--round"), type="character", default="both",
                help="Filter by round (1, 2, or both) [default=%default]"),
    make_option(c("--remove_na"), type="logical", default=TRUE,
                help="Remove NA values from all columns"),
    make_option(c("--remove_empty"), type="logical", default=TRUE,
                help="Remove empty strings from all columns"),
    make_option(c("--remove_null"), type="logical", default=TRUE,
                help="Remove NULL values from all columns"),
    make_option(c("--format"), type="character", default="png",
                help="Output format (png or svg)"),
    make_option(c("--width"), type="integer", default=8000,
                help="Plot width in pixels"),
    make_option(c("--height"), type="integer", default=6000,
                help="Plot height in pixels"),
    make_option(c("--res"), type="integer", default=600,
                help="Plot resolution (DPI)"),
    make_option(c("--main_title"), type="character", default="Correlation of Protein Metrics with Binding Affinity",
                help="Main title for the combined plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Format metric name function 
format_metric_name <- function(metric, for_title = FALSE) {
    # Dictionary of metric name formats - plain text versions for titles
    metric_formats <- list(
        "kd" = if(for_title) "binding affinity" else "K[D]",
        "-log10_kd" = if(for_title) "binding affinity" else "-log[10](K[D])",
        "normalized_kd" = if(for_title) "Normalized KD" else "paste('Normalized ', K[D])",
        "-log10_normalized_kd" = if(for_title) "-log10(KD)" else "-log[10](K[D])",
        "iptm" = "ipTM",
        "ipae" = "iPAE",
        "pae_interaction" = "iPAE",
        "esm_pll" = "ESM PLL",
        "normalized_esm_pll" = "Normalized ESM PLL",
        "plddt" = "pLDDT"
    )
    
    # Convert to lower case for matching
    metric_lower <- tolower(metric)
    
    # Return formatted name if it exists, otherwise clean up the default name
    if (metric_lower %in% names(metric_formats)) {
        formatted <- metric_formats[[metric_lower]]
        return(if(!for_title && !grepl("^-log", metric_lower)) sprintf("'%s'", formatted) else formatted)
    } else {
        # Handle -log10 transformed metrics that aren't in the dictionary
        if (grepl("^-log10_", metric_lower)) {
            base_metric <- sub("^-log10_", "", metric_lower)
            base_format <- format_metric_name(base_metric, for_title)
            # Remove quotes if they exist in base format
            base_format <- gsub("^'|'$", "", base_format)
            if (for_title) {
                return(sprintf("-log10(%s)", base_format))
            } else {
                return(sprintf("-log[10](%s)", base_format))
            }
        }
        # Default cleanup for unknown metrics
        formatted <- stringr::str_to_title(gsub("_", " ", metric))
        return(if(for_title) formatted else sprintf("'%s'", formatted))
    }
}

# Function to create a single correlation plot
create_correlation_plot <- function(data, x_column, y_column, color_by, color_values) {
    # Calculate correlations and p-values
    pearson_test <- suppressWarnings(cor.test(data[[x_column]], data[[y_column]], 
                                            method = "pearson", use = "complete.obs"))
    
    # Format p-values with appropriate notation
    format_p_value <- function(p) {
        if (p < 0.001) return(sprintf("p < 0.001"))
        if (p < 0.01) return(sprintf("p = %.3f", p))
        return(sprintf("p = %.2f", p))
    }
    
    # Create simple correlation text
    cor_text <- sprintf(
        "r = %.3f (%s)",
        pearson_test$estimate, format_p_value(pearson_test$p.value)
    )
    
    # Format axis labels
    x_label <- format_metric_name(x_column, for_title = FALSE)
    y_label <- format_metric_name(y_column, for_title = FALSE)
    x_title <- format_metric_name(x_column, for_title = TRUE)
    
    # Create plot
    p <- ggplot(data, aes(x = get(x_column), 
                         y = get(y_column),
                         color = get(color_by))) +
        geom_point(size = 1.5, alpha = 0.7) +
        suppressWarnings(geom_smooth(method = "lm", se = TRUE, color = "#3E6175", alpha = 0.2)) +
        scale_color_manual(values = color_values) +
        labs(title = x_title,
             subtitle = cor_text,
             x = parse(text = x_label),
             y = parse(text = y_label),
             color = stringr::str_to_sentence(gsub("_", " ", color_by))) +
        adaptyv_theme_no_grid()
    
    return(p)
}

main <- function() {
    tryCatch({
        # Read data
        log_info(sprintf("Reading data from %s", opt$input))
        data <- suppressWarnings(fread(opt$input))
        
        # Filter by round if specified
        if (opt$round != "both") {
            round_num <- as.numeric(opt$round)
            data <- data[round == round_num]
            log_info(sprintf('Filtered data for round %d', round_num))
        }
        
        # Set y_column and transform KD values if needed
        y_column <- opt$y_column
        if (grepl("kd", tolower(y_column))) {
            # Create the new column name with correct prefix
            transformed_y_column <- paste0("-log10_", y_column)
            # Add the transformed column to the data
            data[[transformed_y_column]] <- -log10(data[[y_column]])
            # Update the y_column to use the transformed column
            y_column <- transformed_y_column
            log_info(sprintf("Transformed KD values to %s", y_column))
        }
        
        # Define the x metrics for the 4 subplots
        x_metrics <- c("pae_interaction", "iptm", "esm_pll", "normalized_esm_pll")
        
        # Check if all required columns exist
        required_cols <- c(x_metrics, y_column, opt$color_by)
        if (!all(required_cols %in% names(data))) {
            missing_cols <- setdiff(required_cols, names(data))
            stop(sprintf("Missing columns in data: %s", paste(missing_cols, collapse = ", ")))
        }
        
        # Select appropriate color palette based on color_by column
        color_values <- switch(opt$color_by,
                             "binding_strength" = binding_colors,
                             "expression" = expression_colors,
                             "selection_status" = selection_status_colors,
                             "design_category" = design_category_colors,
                             adaptyv_colors)
        
        # Create a list to store the plots
        plots_list <- list()
        
        # Filter data and create each plot
        for (x_metric in x_metrics) {
            # Create filtering conditions based on options
            filter_conditions <- list()
            
            # Handle NA values
            if (opt$remove_na) {
                filter_conditions[[length(filter_conditions) + 1]] <- 
                    sprintf("!is.na(%s) & !is.na(%s) & !is.na(%s)", 
                           x_metric, y_column, opt$color_by)
            }
            
            # Handle empty strings
            if (opt$remove_empty) {
                filter_conditions[[length(filter_conditions) + 1]] <- 
                    sprintf("%s != '' & %s != '' & %s != ''", 
                           x_metric, y_column, opt$color_by)
            }
            
            # Handle NULL values
            if (opt$remove_null) {
                filter_conditions[[length(filter_conditions) + 1]] <- 
                    sprintf("!is.null(%s) & !is.null(%s) & !is.null(%s)", 
                           x_metric, y_column, opt$color_by)
            }
            
            # Apply filtering
            if (length(filter_conditions) > 0) {
                filter_expr <- paste(filter_conditions, collapse = " & ")
                valid_data <- data[eval(parse(text = filter_expr))]
                log_info(sprintf("Filtered data for %s plot: %d rows", x_metric, nrow(valid_data)))
            } else {
                valid_data <- data
            }
            
            # Create the plot
            plots_list[[x_metric]] <- create_correlation_plot(
                valid_data, x_metric, y_column, opt$color_by, color_values
            )
        }
        
        # Arrange the 4 plots in a 2x2 grid with a common legend
        g <- arrangeGrob(
            plots_list[["pae_interaction"]] + theme(legend.position = "none"),
            plots_list[["iptm"]] + theme(legend.position = "none"),
            plots_list[["esm_pll"]] + theme(legend.position = "none"),
            plots_list[["normalized_esm_pll"]] + theme(legend.position = "none"),
            nrow = 2,
            top = textGrob(opt$main_title, gp = gpar(fontsize = 16, fontface = "bold", fontfamily = "Roboto")),
            bottom = get_legend(plots_list[[1]] + theme(legend.position = "bottom", 
                                                     legend.box = "horizontal",
                                                     legend.margin = margin(t = 10)))
        )
        
        # Save plot
        dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
        filename <- file.path(opt$output, 
                            sprintf("combined_metrics_vs_%s_%s.%s",
                                    gsub("-log10_", "", y_column), opt$color_by, opt$format))
        
        if (opt$format == "png") {
            png(filename, width = opt$width, height = opt$height, res = opt$res)
        } else {
            svg(filename, width = opt$width/150, height = opt$height/150)
        }
        
        grid.draw(g)
        dev.off()
        
        log_info(sprintf("Successfully created combined metrics plot: %s", filename))
        
    }, error = function(e) {
        log_error(sprintf("Error in combined metrics plot creation: %s", e$message))
        stop(e)
    })
}

if (sys.nframe() == 0) {
    main()
} 