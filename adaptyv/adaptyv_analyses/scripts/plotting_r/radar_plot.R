#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(fmsb)  # For radar plots

font_add_google("Roboto", "Roboto")
showtext_auto()

# Reuse the existing theme from other scripts
adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 12, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box.background = element_rect(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Define command line options
option_list <- list(
  make_option(c("--input"), type="character", 
              default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", 
              default="./plots/radar_plots",
              help="Output directory path"),
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
  make_option(c("--binder_type"), type="character", default="both",
              help="Filter by binder type (De novo, Existing binder, or both)"),
  make_option(c("--top_n"), type="integer", default=5,
              help="Number of top binders to plot"),
  make_option(c("--title"), type="character", 
              default="Interface metrics for the top binders",
              help="Custom title for the plot"),
  make_option(c("--subtitle"), type="character", 
              default="Comparing De novo and Existing binders",
              help="Custom subtitle for the plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Interface metrics to plot (now including KD)
interface_metrics <- c(
  "interface_dG",
  "interface_hydrophobicity",
  "interface_nres",
  "interface_hbond_percentage",
  "interface_sc",
  "kd"  # Added KD as a metric
)

# Define color palettes for different binder types with more distinct colors
binder_colors <- list(
  "De novo" = colorRampPalette(c("#56A6D4", "#2A4DD0", "#3E6175", "#8B90DD", "#8CD2F4"))(5),  # Blue spectrum
  "Existing binder" = colorRampPalette(c("#FFB347", "#FF8C00", "#FF6B6B", "#FF4500", "#DC582A"))(5)  # Orange-red spectrum
)

# Function to scale metrics between 0 and 1
scale_metric <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

main <- function() {
  tryCatch({
    # Create output directory if it doesn't exist
    dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
    
    if (!dir.exists(opt$output)) {
      stop(sprintf("Failed to create output directory: %s", opt$output))
    }
    
    # Read and prepare data
    df <- fread(opt$input)
    
    # Filter by round if specified
    if (opt$round != "both") {
      round_num <- as.numeric(opt$round)
      df <- df[round == round_num]
      log_info(sprintf('Filtered data for round %d', round_num))
    }
    
    # Create separate plots for De novo and Existing binders
    binder_types <- if(opt$binder_type == "both") {
      c("De novo", "Existing binder")
    } else {
      c(opt$binder_type)
    }
    
    for(binder_type in binder_types) {
      # Filter by binder type
      df_filtered <- df[de_novo == binder_type]
      
      # Select top N binders by KD
      top_binders <- df_filtered[order(kd)][1:min(opt$top_n, nrow(df_filtered))]
      
      if (nrow(top_binders) == 0) {
        log_warn(sprintf("No data points for %s binders", binder_type))
        next
      }
      
      # Prepare radar plot data
      radar_data <- top_binders[, interface_metrics, with = FALSE]
      radar_data <- as.data.frame(radar_data)
      
      # Create design names with KD values
      design_names <- sprintf("Design %d - KD: %.2f nM", 
                            seq_len(nrow(radar_data)), 
                            top_binders$kd)
      rownames(radar_data) <- design_names
      
      # Scale the metrics (special handling for KD)
      scaled_data <- as.data.frame(lapply(radar_data, function(x) {
        if (identical(x, radar_data$kd)) {
          # For KD, lower is better, so reverse the scaling
          return(1 - scale_metric(x))
        } else {
          return(scale_metric(x))
        }
      }))
      
      # Add max/min rows required by fmsb
      scaled_data <- rbind(rep(1, ncol(scaled_data)),
                          rep(0, ncol(scaled_data)),
                          scaled_data)
      
      # Get appropriate color palette based on binder type
      colors <- binder_colors[[binder_type]]
      
      # Set up plot title and subtitle with proper capitalization
      title <- if(opt$title == "") {
        "Interface metrics for the top binders"  # Only first word capitalized
      } else {
        opt$title
      }
      
      subtitle <- sprintf("%s binders", # Only first word capitalized
                         ifelse(binder_type == "De novo", "De novo", "Existing"))
      
      # Format KD values (already in nM, no conversion needed)
      kd_values <- sapply(top_binders$kd, function(x) {
        formatC(x, format = "e", digits = 1)
      })
      
      # Create design labels with KD values in nM
      design_labels <- sprintf("Design %d (KD: %s nM)", 
                             seq_len(nrow(radar_data)), 
                             kd_values)
      
      # Add metric labels with better formatting
      metric_labels <- c(
        "Interface ΔG",
        "Hydrophobicity",
        "Interface size",
        "H-bond %",
        "Shape comp.",
        "Binding affinity"  # Changed from KD (nM)
      )
      
      # Define min/max values for each metric for scale labels
      metric_ranges <- list(
        "Interface ΔG" = c("Less stable", "More stable"),
        "Hydrophobicity" = c("Hydrophilic", "Hydrophobic"),
        "Interface size" = c("Smaller", "Larger"),
        "H-bond %" = c("Fewer", "More"),
        "Shape comp." = c("Lower", "Higher"),
        "Binding affinity" = c("Lower", "Higher")
      )
      
      # Define line types for each design
      line_types <- c("solid", "dashed", "dotted", "twodash", "longdash")
      
      # Generate filename
      filename <- file.path(opt$output, sprintf(
        "radar_plot_%s_r%s.%s",
        gsub(" ", "_", tolower(binder_type)),
        opt$round,
        opt$format
      ))
      
      # Create plot with error handling
      tryCatch({
        if (opt$format == "png") {
          png(filename, width = opt$width, height = opt$height, res = opt$res)
        } else {
          svg(filename)
        }
        
        par(mar = c(2, 2, 4, 2))  # Adjust margins for title and subtitle
        
        # Create radar chart with clean axis labels and no numerical scale
        radarchart(
          scaled_data,
          pcol = colors,
          plty = line_types[1:nrow(radar_data)],  # Add different line types
          pfcol = NA,
          plwd = 2,
          cglcol = "grey80",
          cglty = 1,
          axislabcol = "grey30",
          caxislabels = rep("", 5),
          axistype = 1,
          seg = 4,
          cex.axis = 0.7,
          cex.lab = 0.8,
          vlabels = metric_labels,
          title = title
        )
        
        # Add subtitle
        mtext(subtitle, 
              side = 3, 
              line = 0.5, 
              cex = 0.8, 
              col = "#3D7A9A")
        
        # Add legend with title, box, and line types
        legend(
          "topright",
          title = expression(bold("Designs ranked by KD")),
          legend = design_labels,
          col = colors,
          lty = line_types[1:length(design_labels)],  # Add line types to legend
          lwd = 2,
          pch = NA,
          bty = "o",
          box.col = "black",
          box.lwd = 0.5,
          cex = 0.7,
          text.col = "grey30"
        )
        
        # Close the device to save the plot
        dev.off()
        log_info(sprintf('Saved plot to %s', filename))
        
        # Save data used for plotting
        data_filename <- sub(paste0(".", opt$format), ".csv", filename)
        write.csv(
          cbind(name = rownames(scaled_data), scaled_data),
          file = data_filename,
          row.names = FALSE
        )
        log_info(sprintf('Saved plot data to %s', data_filename))
        
      }, error = function(e) {
        # Clean up if plot creation fails
        if (dev.cur() > 1) dev.off()
        log_error(sprintf("Failed to create plot: %s", e$message))
        stop(e)
      })
    }
  }, error = function(e) {
    log_error(sprintf("Script failed: %s", e$message))
    stop(e)
  })
}

if (sys.nframe() == 0) {
  main()
}
