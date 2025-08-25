#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr, warn.conflicts = FALSE)
library(logger)
library(optparse)

font_add_google("Roboto", "Roboto")
showtext_auto()

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

adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, 
                                 selection_status_colors, design_category_colors)))

# Add round colors from barplots.R
round_colors <- c(
    "2" = "#8B90DD",  # purple
    "1" = "#8CD2F4"   # blue
)

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.spacing = unit(2, "lines"),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      plot.margin = margin(t = 20, r = 100, b = 20, l = 20, unit = "pt")
    )
}

format_metric_name <- function(metric, for_title = FALSE) {
    metric_formats <- list(
        "kd" = if(for_title) "binding affinity" else "'K'[D]",
        "-log10_kd" = if(for_title) "binding affinity" else "-log[10]('K'[D])",
        "normalized_kd" = if(for_title) "Normalized KD" else "paste('Normalized ', 'K'[D])",
        "-log10_normalized_kd" = if(for_title) "-log10(KD)" else "-log[10]('K'[D])",
        "iptm" = "ipTM",
        "ipae" = "iPAE",
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
        formatted <- stringr::str_to_title(gsub("_", " ", metric))
        return(if(for_title) formatted else sprintf("'%s'", formatted))
    }
}

option_list <- list(
    make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
                help="Input CSV file path"),
    make_option(c("--output"), type="character", default="./plots/violin/",
                help="Output directory path"),
    make_option(c("--y_column"), type="character", default="kd",
                help="Column name for y-axis (metric)"),
    make_option(c("--x_column"), type="character", default="round",
                help="Column name for x-axis grouping"),
    make_option(c("--color_by"), type="character", default="round",
                help="Optional column name for color grouping"),
    make_option(c("--round"), type="character", default='both',
                help="Filter by round (1, 2, or both)"),
    make_option(c("--format"), type="character", default="svg",
                help="Output format (png or svg)"),
    make_option(c("--width"), type="integer", default=1600,
                help="Plot width in pixels"),
    make_option(c("--height"), type="integer", default=1200,
                help="Plot height in pixels"),
    make_option(c("--res"), type="integer", default=300,
                help="Plot resolution (DPI)"),
    make_option(c("--title"), type="character", default="",
                help="Custom title for the plot"),
    make_option(c("--subtitle"), type="character", default="",
                help="Custom subtitle for the plot"),
    make_option(c("--remove_not_mentioned"), type="logical", default=TRUE,
                help="Remove 'Not mentioned' category [default=%default]"),
    make_option(c("--show_anova"), type="logical", default=TRUE,
                help="Show ANOVA test results in subtitle [default=%default]"),
    make_option(c("--binders_only"), type="logical", default=FALSE,
                help="Filter for binders only [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

main <- function() {
    tryCatch({
        log_info('Starting script execution')
        
        log_info(sprintf('Loading data from %s', opt$input))
        data <- fread(opt$input)
        
        if (opt$binders_only) {
            data <- data[!is.na(binding) & binding != "" & binding != "NULL"]
            data <- data[binding == "Yes"]
            log_info('Filtered data for binders only')
        }
        
        if (opt$round != "both") {
            round_num <- as.numeric(opt$round)
            data <- data[round == round_num]
            log_info(sprintf('Filtered data for round %d', round_num))
        }
        
        y_label <- format_metric_name(opt$y_column, for_title = FALSE)
        x_label <- stringr::str_to_title(gsub("_", " ", opt$x_column))
        
        if (grepl("kd", tolower(opt$y_column))) {
            data[[paste0("-log10_", opt$y_column)]] <- -log10(data[[opt$y_column]])
            opt$y_column <- paste0("-log10_", opt$y_column)
            y_label <- format_metric_name(opt$y_column, for_title = FALSE)
        }
        
        valid_data <- data[!is.na(get(opt$y_column)) & !is.na(get(opt$x_column))]
        if (!is.null(opt$color_by)) {
            valid_data <- valid_data[!is.na(get(opt$color_by))]
        }
        
        if (opt$remove_not_mentioned) {
            valid_data <- valid_data[get(opt$x_column) != "Not mentioned"]
            if (!is.null(opt$color_by)) {
                valid_data <- valid_data[get(opt$color_by) != "Not mentioned"]
            }
            log_info('Removed "Not mentioned" category')
        }
        if (x_label == "selected") {
            x_label <- "selection status"
        }
        if (opt$title == "") {
            title <- sprintf("Distribution of %s by %s",
                           format_metric_name(opt$y_column, for_title = TRUE),
                           tolower(x_label))
        } else {
            title <- opt$title
        }
        
        subtitle <- ""
        if (opt$show_anova) {
            # Replace ANOVA with Kruskal-Wallis test
            kw_result <- kruskal.test(get(opt$y_column) ~ get(opt$x_column), data = valid_data)
            
            # Calculate effect size (epsilon-squared)
            n <- nrow(valid_data)
            eps_squared <- (kw_result$statistic - length(unique(valid_data[[opt$x_column]])) + 1) / (n - 1)
            
            # Format p-value
            if (kw_result$p.value < 0.001) {
                p_formatted <- "p < 0.001"
            } else {
                p_formatted <- sprintf("p = %.3f", kw_result$p.value)
            }
            
            subtitle <- sprintf("Kruskal-Wallis: H(%d) = %.2f, %s",
                              kw_result$parameter,
                              kw_result$statistic,
                              p_formatted)
            
            # Add pairwise Wilcoxon tests if significant
            if (kw_result$p.value < 0.05) {
                pairwise_tests <- pairwise.wilcox.test(
                    valid_data[[opt$y_column]], 
                    valid_data[[opt$x_column]],
                    p.adjust.method = "bonferroni"
                )
                # Store results for plotting significance brackets later
                pw_results <- pairwise_tests$p.value
            }
        }
        
        if (opt$subtitle != "") {
            subtitle <- opt$subtitle
        }
        
        #KD control measurements, with data from https://foundry.adaptyvbio.com/competition
        cetuximab_kd <- mean(c(9.94e-9, 3.33e-9))
        egf_kd <- mean(c(7.95e-7, 7.57e-7, 7.25e-7))
        
        # Convert round to factor before plotting
        if (opt$x_column == "round") {
            valid_data$round <- factor(valid_data$round)
            p <- ggplot(valid_data, aes(x = round, y = get(opt$y_column))) +
                geom_violin(aes(fill = round), trim = FALSE, alpha = 0.8) +
                geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA, color = "black", fill = "white") +
                {if(grepl("kd", tolower(opt$y_column))) 
                    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3)
                } +
                stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
                scale_fill_manual(values = round_colors, guide = "none") +
                labs(title = title,
                     subtitle = subtitle,
                     x = x_label,
                     y = parse(text = y_label)) +
                adaptyv_theme() +
                coord_cartesian(clip = "off")
        } else {
            p <- ggplot(valid_data, aes(x = get(opt$x_column), y = get(opt$y_column))) +
                geom_violin(trim = FALSE, alpha = 0.8) +
                geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA, color = "black", fill = "white") +
                {if(grepl("kd", tolower(opt$y_column))) 
                    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3)
                } +
                stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
                labs(title = title,
                     subtitle = subtitle,
                     x = x_label,
                     y = parse(text = y_label)) +
                adaptyv_theme() +
                coord_cartesian(clip = "off")
        }
        
        if (grepl("kd", tolower(opt$y_column))) {
            p <- p +
                geom_hline(yintercept = -log10(cetuximab_kd), 
                          linetype = "dashed", 
                          color = "#E9C435",
                          size = 0.5) +
                geom_hline(yintercept = -log10(egf_kd), 
                          linetype = "dashed", 
                          color = "#69B7A7",
                          size = 0.5) +
                annotate("text", 
                        x = length(unique(valid_data[[opt$x_column]])) + 0.35, 
                        y = -log10(cetuximab_kd) + 0.15, 
                        label = "Cetuximab control", 
                        hjust = 0,  
                        vjust = 0,  
                        size = 3,
                        color = "#E9C435") +
                annotate("text", 
                        x = length(unique(valid_data[[opt$x_column]])) + 0.35, 
                        y = -log10(egf_kd) - 0.15, 
                        label = "EGF control", 
                        hjust = 0,  
                        vjust = 1,  
                        size = 3,
                        color = "#69B7A7") +
                # Add extra space on the right for labels
                scale_x_discrete(expand = expansion(mult = c(0.05, 0.6))) +
                theme(plot.margin = margin(t = 20, r = 120, b = 20, l = 20, unit = "pt")) +
                coord_cartesian(clip = "off")
        }
        
        if (!is.null(opt$color_by)) {
            color_values <- switch(opt$color_by,
                                 "binding_strength" = binding_colors,
                                 "expression" = expression_colors,
                                 "selected" = selection_status_colors,
                                 "design_category" = design_category_colors,
                                 "round" = round_colors,
                                 adaptyv_colors)
            
            p <- p + aes(fill = get(opt$color_by)) +
                scale_fill_manual(values = color_values,
                                name = stringr::str_to_title(gsub("_", " ", opt$color_by)))
        } else {
            p <- p + geom_violin(fill = "#8CD2F4", trim = FALSE, alpha = 0.8) +
                    geom_boxplot(width = 0.1, fill = "white")
        }
        
        dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
        filename <- file.path(opt$output, 
                            sprintf("violin_%s_by_%s%s.%s",
                                    opt$y_column, 
                                    opt$x_column,
                                    ifelse(is.null(opt$color_by), "", sprintf("_%s", opt$color_by)),
                                    opt$format))
        
        if (opt$format == "png") {
            png(filename, width = opt$width, height = opt$height, res = opt$res)
        } else {
            svg(filename, width = opt$width / opt$res, height = opt$height / opt$res)
        }
        print(p)
        dev.off()
        
        log_info(sprintf('Saved plot to %s', filename))
        log_info('Script completed successfully')
        
    }, error = function(e) {
        log_error(sprintf("Script failed with error: %s", e$message))
        stop(e)
    })
}

if (sys.nframe() == 0) {
    main()
} 