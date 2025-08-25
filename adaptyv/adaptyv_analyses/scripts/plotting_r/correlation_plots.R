#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger, warn.conflicts = FALSE)
library(optparse)

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

# Reuse the theme from adaptyv_theme.R
adaptyv_theme <- function() {
    theme_minimal() +
        theme(
            text = element_text(family = "Roboto", color = "#333333"),
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
            plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
            axis.title = element_text(size = 10, face = "bold"),
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
            legend.box.background = element_rect(color = "white"),
            plot.margin = margin(20, 20, 20, 20)
        )
}

# Define command line options
option_list <- list(
    make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
                help="Input CSV file path"),
    make_option(c("--output"), type="character", default="./plots/correlations/",
                help="Output directory path"),
    make_option(c("--x_column"), type="character", default="iptm",
                help="Column name for x-axis"),
    make_option(c("--y_column"), type="character", default="kd",
                help="Column name for y-axis"),
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
    make_option(c("--format"), type="character", default="svg",
                help="Output format (png or svg)"),
    make_option(c("--width"), type="integer", default=4000,
                help="Plot width in pixels"),
    make_option(c("--height"), type="integer", default=3600,
                help="Plot height in pixels"),
    make_option(c("--res"), type="integer", default=600,
                help="Plot resolution (DPI)"),
    make_option(c("--title"), type="character", default=NULL,
                help="Custom plot title"),
    make_option(c("--subtitle"), type="character", default=NULL,
                help="Custom plot subtitle")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Add this function near the top of the script after the color definitions
format_metric_name <- function(metric, for_title = FALSE) {
    # Dictionary of metric name formats - plain text versions for titles
    metric_formats <- list(
        "kd" = if(for_title) "binding affinity" else "K[D]",
        "-log10_kd" = if(for_title) "binding affinity" else "-log[10](K[D])",
        "normalized_kd" = if(for_title) "Normalized KD" else "paste('Normalized ', K[D])",
        "-log10_normalized_kd" = if(for_title) "-log10(KD)" else "-log[10](K[D])",
        "iptm" = "ipTM",
        "ipae" = "iPAE",
        "esm_pll" = "ESM PLL",
        "normalized_esm_pll" = "Normalized ESM PLL",
        "plddt" = "pLDDT",
        "pae_interaction" = "iPAE"
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
        
        # Create filtering conditions based on options
        filter_conditions <- list()
        
        # Handle NA values
        if (opt$remove_na) {
            filter_conditions[[length(filter_conditions) + 1]] <- 
                sprintf("!is.na(%s) & !is.na(%s) & !is.na(%s)", 
                       opt$x_column, opt$y_column, opt$color_by)
            log_info("Removing NA values")
        }
        
        # Handle empty strings
        if (opt$remove_empty) {
            filter_conditions[[length(filter_conditions) + 1]] <- 
                sprintf("%s != '' & %s != '' & %s != ''", 
                       opt$x_column, opt$y_column, opt$color_by)
            log_info("Removing empty strings")
        }
        
        # Handle NULL values
        if (opt$remove_null) {
            filter_conditions[[length(filter_conditions) + 1]] <- 
                sprintf("!is.null(%s) & !is.null(%s) & !is.null(%s)", 
                       opt$x_column, opt$y_column, opt$color_by)
            log_info("Removing NULL values")
        }
        
        # Apply all filtering conditions
        if (length(filter_conditions) > 0) {
            filter_expr <- paste(filter_conditions, collapse = " & ")
            valid_data <- data[eval(parse(text = filter_expr))]
            log_info(sprintf("Filtered from %d to %d rows", nrow(data), nrow(valid_data)))
        } else {
            valid_data <- data
            log_info("No filtering applied")
        }
        
        # Check if columns exist
        if (!all(c(opt$x_column, opt$y_column, opt$color_by) %in% names(valid_data))) {
            missing_cols <- setdiff(c(opt$x_column, opt$y_column, opt$color_by), names(valid_data))
            stop(sprintf("Missing columns in data: %s", paste(missing_cols, collapse = ", ")))
        }

        # Transform KD if necessary BEFORE calculating correlations
        if (grepl("kd", tolower(opt$x_column))) {
            valid_data[[paste0("-log10_", opt$x_column)]] <- -log10(valid_data[[opt$x_column]])
            opt$x_column <- paste0("-log10_", opt$x_column)
            log_info(sprintf("Transformed x-axis column to %s", opt$x_column))
        }
        if (grepl("kd", tolower(opt$y_column))) {
            valid_data[[paste0("-log10_", opt$y_column)]] <- -log10(valid_data[[opt$y_column]])
            opt$y_column <- paste0("-log10_", opt$y_column)
            log_info(sprintf("Transformed y-axis column to %s", opt$y_column))
        }
      
        log_info(sprintf("Removed x-axis outliers: %d rows remaining", nrow(valid_data)))
        
        # Calculate correlations and p-values using potentially transformed columns
        pearson_test <- suppressWarnings(cor.test(valid_data[[opt$x_column]], valid_data[[opt$y_column]], 
                                method = "pearson", use = "complete.obs"))
        spearman_test <- suppressWarnings(cor.test(valid_data[[opt$x_column]], valid_data[[opt$y_column]], 
                                 method = "spearman", use = "complete.obs"))
        kendall_test <- suppressWarnings(cor.test(valid_data[[opt$x_column]], valid_data[[opt$y_column]], 
                                method = "kendall", use = "complete.obs"))
        
        # Format p-values with appropriate notation
        format_p_value <- function(p) {
            if (p < 0.001) return(sprintf("p < 0.001"))
            if (p < 0.01) return(sprintf("p = %.3f", p))
            return(sprintf("p = %.2f", p))
        }
        
        # Create correlation text with colored statistics
        cor_text <- sprintf(
            "<b>Correlation coefficients:</b><br>Pearson <i>r</i> = <span style='color:#3D7A9A'>%.3f</span> (%s)<br>Spearman ρ = <span style='color:#3D7A9A'>%.3f</span> (%s)<br>Kendall τ = <span style='color:#3D7A9A'>%.3f</span> (%s)",
            pearson_test$estimate, format_p_value(pearson_test$p.value),
            spearman_test$estimate, format_p_value(spearman_test$p.value),
            kendall_test$estimate, format_p_value(kendall_test$p.value)
        )
        
        # Set plot title
        if (is.null(opt$title)) {
            title <- sprintf("%s vs %s",
                           gsub("_", " ", opt$x_column),
                           gsub("_", " ", opt$y_column))
        } else {
            title <- opt$title
        }
        
        # Use correlation text as subtitle if no custom subtitle provided
        if (is.null(opt$subtitle)) {
            subtitle <- cor_text
        } else {
            subtitle <- opt$subtitle
        }
        
        # Select appropriate color palette based on color_by column
        color_values <- switch(opt$color_by,
                             "binding_strength" = binding_colors,
                             "expression" = expression_colors,
                             "selection_status" = selection_status_colors,
                             "design_category" = design_category_colors,
                             adaptyv_colors)
        
        # Create base plot
        # KD transformation is now done before correlation calculation, so remove from here
        # if (grepl("kd", tolower(opt$x_column))) {
        #     valid_data[[paste0("-log10_", opt$x_column)]] <- -log10(valid_data[[opt$x_column]])
        #     opt$x_column <- paste0("-log10_", opt$x_column)
        # }
        # if (grepl("kd", tolower(opt$y_column))) {
        #     valid_data[[paste0("-log10_", opt$y_column)]] <- -log10(valid_data[[opt$y_column]])
        #     opt$y_column <- paste0("-log10_", opt$y_column)
        # }

        # Format axis labels and title differently
        x_label <- format_metric_name(opt$x_column, for_title = FALSE)
        y_label <- format_metric_name(opt$y_column, for_title = FALSE)
        x_title <- format_metric_name(opt$x_column, for_title = TRUE)
        y_title <- format_metric_name(opt$y_column, for_title = TRUE)
        
        # Update title with formatted metric names
        if (is.null(opt$title)) {
            title <- sprintf("%s vs %s",
                           x_title, y_title)
        } else {
            title <- opt$title
        }
        
        # Create plot with formatted labels
        p <- ggplot(valid_data, aes(x = get(opt$x_column), 
                                   y = get(opt$y_column),
                                   color = get(opt$color_by))) +
            geom_point(size = 2) +
            suppressWarnings(geom_smooth(method = "lm", se = TRUE, color = "#3E6175", alpha = 0.2)) +
            scale_color_manual(values = color_values) +
            labs(title = title,
                 subtitle = subtitle,
                 x = parse(text = x_label),
                 y = parse(text = y_label),
                 color = stringr::str_to_sentence(gsub("_", " ", opt$color_by))) +
            adaptyv_theme()
            
       
        # Save plot
        dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
        filename <- file.path(opt$output, 
                            sprintf("correlation_%s_vs_%s_%s.%s",
                                    opt$x_column, opt$y_column, opt$color_by, opt$format))
        
        if (opt$format == "png") {
            png(filename, width = opt$width, height = opt$height, res = opt$res)
        } else {
            svg(filename, width = opt$width/opt$res, height = opt$height/opt$res)
        }
        print(p)
        dev.off()
        
        log_info(sprintf("Successfully created correlation plot: %s", filename))
        
    }, error = function(e) {
        log_error(sprintf("Error in correlation plot creation: %s", e$message))
        stop(e)
    })
}

if (sys.nframe() == 0) {
    main()
}