#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(grid)
library(gridExtra)

font_add_google("Roboto", "Roboto")
showtext_auto()

# Interface metrics to plot (in desired order)
interface_metrics <- c(
  "interface_hbond_percentage",
  "interface_delta_unsat_hbonds_percentage",
  "interface_dG",
  "interface_nres",
  "interface_hydrophobicity",
  "interface_sc"
)

# Define color palette for de novo vs existing binders vs non-binders
binding_status_colors <- c(
  "De novo" = "#67AFD8",
  "Existing binder" = "#8D92DE",
  "Non-binder" = "#E5E7EB"
)

# Command line options
option_list <- list(
  make_option(c("--input"), type="character", 
              default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", 
              default="./plots/interface_violins/",
              help="Output directory path"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=5000,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=3800,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI)"),
  make_option(c("--round"), type="character", default="both",
              help="Filter by round (1, 2, or both)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

remove_outliers <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm=TRUE)
  H <- 10 * IQR(x, na.rm=TRUE)  
  x[x < (qnt[1] - H) | x > (qnt[2] + H)] <- NA
  return(x)
}

# Function to calculate statistics and create label
calculate_stats <- function(data, metric_name) {
  # Filter data for this metric
  metric_data <- data[data$metric == metric_name, ]
  metric_data <- metric_data[!is.na(metric_data$value), ]
  
  # Run Kruskal-Wallis test
  kw_test <- kruskal.test(value ~ design_category, data = metric_data)
  
  # Get statistics
  h_stat <- kw_test$statistic
  p_val <- kw_test$p.value
  df <- kw_test$parameter
  
  # Format p-value text
  p_text <- if(p_val < 0.001) "p < 0.001"
            else sprintf("p = %.3f", p_val)
  
  # Get y-axis range
  y_max <- max(metric_data$value, na.rm = TRUE)
  y_min <- min(metric_data$value, na.rm = TRUE)
  y_range <- y_max - y_min
  
  # Run pairwise Wilcoxon test with Bonferroni correction
  pairwise_test <- pairwise.wilcox.test(metric_data$value, 
                                       metric_data$design_category,
                                       p.adjust.method = "bonferroni")
  
  # Convert matrix to data frame
  p_values <- as.data.frame(as.table(pairwise_test$p.value))
  names(p_values) <- c("group1", "group2", "p.adj")
  p_values$comparison <- paste(p_values$group1, p_values$group2, sep="-")
  
  list(
    anova_stats = data.frame(
      metric = metric_name,
      x = 2,
      y = y_max + 0.15 * y_range,
      label = sprintf("H(%d) = %.2f, %s", df, h_stat, p_text)
    ),
    tukey = p_values,
    y_max = y_max,
    y_range = y_range
  )
}

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
    
    # Create design_category column
    df[, design_category := ifelse(binding == "No", "Non-binder",
                                  ifelse(binding == "Yes" & de_novo == "De novo", "De novo",
                                        ifelse(binding == "Yes" & de_novo == "Existing binder", "Existing binder", NA)))]
    
    # Filter out NA and Unknown binding status
    df <- df[!is.na(design_category)]
    
    # Remove outliers and prepare plot data
    plot_data <- data.table()
    for (metric in interface_metrics) {
      # Remove outliers
      df[[metric]] <- remove_outliers(df[[metric]])
      
      # Convert design_category to factor with specific order
      df$design_category <- factor(df$design_category, 
                                 levels = c("De novo", "Existing binder", "Non-binder"))
      
      # Prepare long format data
      metric_data <- data.table(
        design_category = df$design_category,
        metric = metric,
        value = df[[metric]]
      )
      plot_data <- rbind(plot_data, metric_data)
    }
    
    # Create better metric labels
    metric_labels <- c(
      "interface_hbond_percentage" = "H-bond %",
      "interface_delta_unsat_hbonds_percentage" = "Δ Unsat. H-bonds %",
      "interface_dG" = "Interface ΔG",
      "interface_nres" = "Interface size",
      "interface_hydrophobicity" = "Hydrophobicity",
      "interface_sc" = "Shape comp."
    )
    
    # Convert metric to factor with specific order
    plot_data$metric <- factor(plot_data$metric, 
                             levels = c("interface_hbond_percentage",
                                      "interface_delta_unsat_hbonds_percentage",
                                      "interface_dG",
                                      "interface_nres",
                                      "interface_hydrophobicity",
                                      "interface_sc"))
    
    # Calculate statistics for all metrics
    all_stats <- lapply(interface_metrics, calculate_stats, data = plot_data)
    names(all_stats) <- interface_metrics

    # Create base plot
    p <- ggplot(plot_data, aes(x = design_category, y = value, fill = design_category)) +
      # Add jittered points first (so they appear behind)
      geom_jitter(data = plot_data, # Only show points for binders
                 width = 0.2, size = 0.3, alpha = 0.2) +
      geom_violin(trim = FALSE, alpha = 0.8) +
      # Add median and quartiles with white fill
      geom_boxplot(width = 0.2, alpha = 0.7, outlier.shape = NA, color = "black", fill = "white") +
      # Add mean point - on top and matching subtitle blue
      stat_summary(fun = mean, geom = "point", size = 1, color = "black") +
      facet_wrap(~metric, 
                labeller = labeller(metric = metric_labels),
                scales = "free_y",
                nrow = 2) +
      scale_fill_manual(values = binding_status_colors, name = "Binder type: ") +
      labs(title = "Interface properties comparison",
           subtitle = "Distribution of values across binder categories",
           x = "Design category",
           y = "Value")

    # Add Kruskal-Wallis statistics
    anova_stats <- do.call(rbind, lapply(all_stats, function(x) x$anova_stats))
    p <- p + geom_text(
      data = anova_stats,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 2.5,
      color = "#3D7A9A",
      fontface = "bold"
    )

    # Add significance stars for all metrics
    for(m in names(all_stats)) {
      stats <- all_stats[[m]]
      tukey_results <- stats$tukey
      
      # Process each pairwise comparison
      for(i in 1:nrow(tukey_results)) {
        p.adj <- tukey_results$p.adj[i]
        if(!is.na(p.adj) && p.adj < 0.05) {  # Check for NA before comparison
          # Get comparison groups
          comp <- tukey_results$comparison[i]
          groups <- strsplit(comp, "-")[[1]]
          x1 <- which(levels(plot_data$design_category) == groups[1])
          x2 <- which(levels(plot_data$design_category) == groups[2])
          
          # Determine significance level
          stars <- if(p.adj < 0.001) "***"
                  else if(p.adj < 0.01) "**"
                  else "*"
          
          # Add significance stars
          p <- p + geom_text(
            data = data.frame(
              metric = m,
              x = mean(c(x1, x2)),
              y = stats$y_max - 0.08 * stats$y_range  # Adjusted position
            ),
            aes(x = x, y = y),
            label = stars,
            inherit.aes = FALSE,
            color = "#3D7A9A",
            size = 3
          )
        }
      }
    }

    # Add theme
    p <- p + theme_minimal() +
      theme(
        text = element_text(family = "Roboto", color = "#333333"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "right",
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid = element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
      ) +
      coord_cartesian(clip = "off")
    
    # Save plot
    dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
    filename <- file.path(opt$output, sprintf("interface_violins_round%s.%s", opt$round, opt$format))
    
    if (opt$format == "png") {
      png(filename, width = opt$width, height = opt$height, res = opt$res)
    } else {
      svg(filename, width = opt$width / opt$res, height = opt$height / opt$res)
    }
    
    print(p)
    dev.off()
    
    log_info(sprintf("Plot saved to %s", filename))
    
  }, error = function(e) {
    log_error(sprintf("Script failed with error: %s", e$message))
    stop(e)
  })
}

if (sys.nframe() == 0) {
  main()
}
