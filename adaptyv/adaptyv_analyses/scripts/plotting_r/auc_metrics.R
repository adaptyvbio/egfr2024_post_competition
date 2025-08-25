#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(pROC)

font_add_google("Roboto", "Roboto")
showtext_auto()

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.box.background = element_rect(color = "black"),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(20, 20, 20, 20)
    )
}

option_list <- list(
  make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", default="./plots/roc_curves",
              help="Output directory path"),
  make_option(c("--metric"), type="character", default="esm_pll",
              help="Metric to analyze (e.g., iptm, pae_interaction, esm_pll)"),
  make_option(c("--type"), type="character", default="both",
              help="Type of ROC curve: 'binding', 'expression', or 'both'"),
  make_option(c("--format"), type="character", default="png",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=3600,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=3600,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI)"),
  make_option(c("--round"), type="character", default=2,
              help="Filter by round (1, 2, or both)"),
  make_option(c("--title"), type="character", default=NULL,
              help="Custom plot title"),
  make_option(c("--subtitle"), type="character", default=NULL,
              help="Custom plot subtitle")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

format_metric_name <- function(metric, for_title = FALSE) {
    metric_formats <- list(
        "kd" = if(for_title) "KD" else "'K'[D]",
        "-log10_kd" = if(for_title) "-log10(KD)" else "-log[10]('K'[D])",
        "normalized_kd" = if(for_title) "Normalized KD" else "paste('Normalized ', 'K'[D])",
        "-log10_normalized_kd" = if(for_title) "-log10(KD)" else "-log[10]('K'[D])",
        "iptm" = "ipTM",
        "pae_interaction" = "iPAE",
        "esm_pll" = "ESM PLL",
        "normalized_esm_pll" = "Normalized ESM PLL",
        "plddt" = "pLDDT",
        "pae_score" = "PAE Score"
    )
    
    metric_lower <- tolower(metric)
    if (metric_lower %in% names(metric_formats)) {
        formatted <- metric_formats[[metric_lower]]
        return(if(!for_title && !grepl("^-log", metric_lower)) sprintf("'%s'", formatted) else formatted)
    } else {
        if (grepl("^-log10_", metric_lower)) {
            base_metric <- sub("^-log10_", "", metric_lower)
            base_format <- format_metric_name(base_metric, for_title)
            base_format <- gsub("^'|'$", "", base_format)
            return(if(for_title) sprintf("-log10(%s)", base_format) 
                   else sprintf("-log[10](%s)", base_format))
        }
        formatted <- gsub("_", " ", metric)
        return(if(for_title) formatted else sprintf("'%s'", formatted))
    }
}

main <- function() {
  df <- fread(opt$input)
  
  if (opt$round != "both") {
    tryCatch({
      round_num <- as.numeric(opt$round)
      if (is.na(round_num)) {
        stop("Invalid round number provided")
      }
      df <- df[round == round_num]
      log_info(sprintf('Filtered data for round %d', round_num))
    }, error = function(e) {
      log_error(sprintf("Error filtering by round: %s", e$message))
      stop(e)
    })
  }
  
  calculate_roc <- function(data, success_col, metric_col) {
    tryCatch({
      valid_data <- data[!is.na(get(metric_col)) & !is.na(get(success_col))]
      if (nrow(valid_data) == 0) {
        stop("No valid data points after filtering NAs")
      }
      
      # Print summary statistics for debugging
      log_info(sprintf("Number of valid data points: %d", nrow(valid_data)))
      log_info(sprintf("Success column unique values: %s", 
                      paste(unique(valid_data[[success_col]]), collapse=", ")))
      log_info(sprintf("Metric range: [%.2f, %.2f]", 
                      min(valid_data[[metric_col]]), max(valid_data[[metric_col]])))
      
      # Set direction based on metric type
      direction <- if(grepl("^iptm$|^plddt$", metric_col, ignore.case=TRUE)) "<" 
                  else if(grepl("^pae_interaction$|^sequence_length$", metric_col, ignore.case=TRUE)) ">"
                  else if(grepl("^esm_pll$|^normalized_esm_pll$", metric_col, ignore.case=TRUE)) "auto" # "<"
                  else "auto"

      roc_obj <- roc(valid_data[[success_col]], valid_data[[metric_col]], 
                     smooth = FALSE, direction = direction)
      coords <- coords(roc_obj, "all")
      data.table(
        specificity = coords$specificity,
        sensitivity = coords$sensitivity,
        auc = auc(roc_obj),
        type = success_col
      )
    }, error = function(e) {
      log_error(sprintf("Error calculating ROC: %s", e$message))
      stop(e)
    })
  }
  
  plot_data <- data.table()
  
  if (opt$type %in% c("binding", "both")) {
    binding_data <- copy(df[!is.na(binding) & binding != "" & binding != "NULL" & binding != "Unknown"]) # Remove NA, empty, NULL, and Unknown (not expressed) values
    binding_data[, binding_success := ifelse(binding == "Yes", 1, 0)]
    plot_data <- rbind(plot_data, calculate_roc(binding_data, "binding_success", opt$metric))
  }
  
  if (opt$type %in% c("expression", "both")) {
    expression_data <- copy(df[!is.na(expression) & expression != "" & expression != "NULL"])
    expression_data[, expression_success := ifelse(expression == "High", 1, 0)]
    plot_data <- rbind(plot_data, calculate_roc(expression_data, "expression_success", opt$metric))
  }

  colors <- c("binding_success" = "#8CD2F4", "expression_success" = "#8B90DD")
  
  title <- ifelse(is.null(opt$title),
                 sprintf("ROC curve for %s", format_metric_name(opt$metric, for_title = TRUE)),
                 opt$title)
  
  subtitle <- ifelse(is.null(opt$subtitle),
                    "Discriminating binding and expression outcomes",
                    opt$subtitle)
  
  p <- ggplot(plot_data, aes(x = 1 - specificity, y = sensitivity, color = type)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = colors,
                      labels = c("binding_success" = "Binding",
                               "expression_success" = "Expression")) +
    geom_text(data = plot_data[, .(x = 0.75, y = ifelse(type == "binding_success", 0.25, 0.15),
                                  label = sprintf("AUC = %.2f", auc),
                                  type = type)],
              aes(x = x, y = y, label = label, color = type),
              size = 5, show.legend = FALSE) +
    labs(x = "False positive rate",
         y = "True positive rate",
         title = title,
         subtitle = subtitle,
         color = "Prediction type",
         size = 10) +
    adaptyv_theme() +
    theme(legend.position = "bottom")
  
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  
  filename <- file.path(opt$output, sprintf("roc_curves_%s_%s.%s", opt$type, opt$metric, opt$format))
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