#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(data.table)
library(showtext)
library(ggtext)
library(dplyr)
library(gridExtra)
library(grid)

# Import theme and colors
font_add_google("Roboto", "Roboto")
showtext_auto()

# Define color palettes
design_category_colors <- c(
  "De novo" = "#56A6D4",
  "Optimized binder" = "#8B90DD",
  "Diversified binder" = "#8CD2F4",
  "Hallucination" = "#3E6175",
  "Not mentioned" = "#F9F5C9"
)

# Helper function to extract legend from a ggplot
get_legend <- function(a_ggplot) {
  tmp <- ggplot_gtable(ggplot_build(a_ggplot + 
                                  theme(legend.position = "right")))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Modified theme without gray grid
theme_no_grid <- function() {
    theme_minimal() +
        theme(
            text = element_text(family = "Roboto", color = "#333333"),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 5)),
            plot.subtitle = element_markdown(size = 6, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 5)),
            axis.title = element_text(size = 9, face = "bold"),
            axis.title.y = element_text(angle = 90, vjust = 0.5),
            axis.text = element_text(size = 7),
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
            legend.position = "right",
            legend.box.background = element_rect(color = "white"),
            plot.margin = margin(3, 3, 3, 3)
        )
}

# Main function
main <- function() {
    # Read data
    cat("Reading data...\n")
    data_file <- "./data/processed/all_submissions.csv"
    data <- fread(data_file)
    
    # Create the log-transformed KD column
    cat("Transforming KD values...\n")
    data$neg_log10_kd <- -log10(data$kd)
    
    # Define the metrics to plot
    metrics <- c("pae_interaction", "iptm", "esm_pll", "normalized_esm_pll")
    metric_names <- c("iPAE", "ipTM", "ESM PLL", "Normalized ESM PLL")
    color_by <- "design_category"
    
    # Create a list to store the plots
    plots_list <- list()
    
    # Create each individual plot
    for (i in 1:length(metrics)) {
        # Filter data for this metric
        valid_data <- data[!is.na(get(metrics[i])) & 
                         !is.na(neg_log10_kd) & 
                         !is.na(get(color_by)) &
                         get(metrics[i]) != "" & 
                         get(color_by) != ""]
        
        cat(sprintf("Creating plot for %s with %d data points\n", 
                   metric_names[i], nrow(valid_data)))
        
        # Calculate correlations
        pearson_result <- cor.test(valid_data[[metrics[i]]], 
                                valid_data$neg_log10_kd, 
                                method = "pearson")
                                
        spearman_result <- cor.test(valid_data[[metrics[i]]], 
                                 valid_data$neg_log10_kd, 
                                 method = "spearman")
                                 
        kendall_result <- cor.test(valid_data[[metrics[i]]], 
                                valid_data$neg_log10_kd, 
                                method = "kendall")
        
        # Format p-values
        format_p_value <- function(p) {
            if (p < 0.001) return("p < 0.001")
            if (p < 0.01) return(sprintf("p = %.3f", p))
            return(sprintf("p = %.2f", p))
        }
        
        # Create correlation text with all three coefficients
        corr_text <- sprintf(
            "Pearson r = %.3f (%s)<br>Spearman ρ = %.3f (%s)<br>Kendall τ = %.3f (%s)",
            pearson_result$estimate, format_p_value(pearson_result$p.value),
            spearman_result$estimate, format_p_value(spearman_result$p.value),
            kendall_result$estimate, format_p_value(kendall_result$p.value)
        )
        
        # Create the plot
        p <- ggplot(valid_data, 
                   aes(x = get(metrics[i]), 
                      y = neg_log10_kd,
                      color = get(color_by))) +
            geom_point(size = 1.2, alpha = 0.7) +
            geom_smooth(method = "lm", se = TRUE, color = "#3E6175", alpha = 0.2) +
            scale_color_manual(values = design_category_colors) +
            labs(title = metric_names[i],
                 subtitle = corr_text,
                 x = metric_names[i],
                 y = "-log10(KD)",
                 color = "Design category") +
            theme_no_grid()
        
        plots_list[[i]] <- p
    }
    
    # Arrange the plots in a 2x2 grid with legend on the right
    cat("Arranging plots...\n")
    g <- grid.arrange(
        arrangeGrob(
            plots_list[[1]] + theme(legend.position = "none"),
            plots_list[[2]] + theme(legend.position = "none"),
            plots_list[[3]] + theme(legend.position = "none"),
            plots_list[[4]] + theme(legend.position = "none"),
            nrow = 2,
            top = textGrob("Correlation of protein metrics with binding affinity",
                         gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Roboto"))
        ),
        get_legend(plots_list[[1]] + theme(legend.position = "right")),
        ncol = 2,
        widths = c(5, 1)
    )
    
    # Save the plot
    output_dir <- "./plots/combined_metrics/"
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(output_dir, "combined_metrics_plot.png")
    
    cat(sprintf("Saving plot to %s\n", output_file))
    png(output_file, width = 6000, height = 4500, res = 600)
    grid.draw(g)
    dev.off()
    
    cat("Done!\n")
}

# Run the main function
main() 