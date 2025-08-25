#!/usr/bin/env Rscript

# Basic print for debugging
print("Script starting...")

suppressPackageStartupMessages({
  library(ggplot2)
  library(showtext)
  library(scales)
  library(ggtext)
  library(dplyr)
  library(optparse)
  library(logger)
})

# Force logger to show messages
log_threshold(INFO)
print("Libraries loaded...")

# Add Roboto font
font_add_google("Roboto", "Roboto")
showtext_auto()

# Theme definition
theme_protein_ml <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 24, hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 16, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 12),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "#E1E8ED", size = 0.5),
      panel.grid.minor = element_line(color = "#E1E8ED", size = 0.25),
      axis.line = element_line(color = "black", size = 0.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(20, 20, 20, 20)
    )
}

# Command line options
option_list <- list(
  make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
              help="Input CSV file path [default=%default]"),
  make_option(c("--output"), type="character", default="./plots/barplots/",
              help="Output directory path [default=%default]"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg) [default=%default]"),
  make_option(c("--width"), type="integer", default=4800,
              help="Plot width in pixels [default=%default]"),
  make_option(c("--height"), type="integer", default=3200,
              help="Plot height in pixels [default=%default]"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI) [default=%default]"),
  make_option(c("--x_column"), type="character", default="model_category",
              help="Column name to plot on x-axis [default=%default]"),
  make_option(c("--color_column"), type="character", default="round",
              help="Column name to use for stacked bars coloring [default=%default]"),
  make_option(c("--remove_not_mentioned"), type="logical", default=FALSE,
              help="Remove 'Not mentioned' category [default=%default]"),
  make_option(c("--remove_missing_data"), type="logical", default=FALSE,
              help="Remove missing data [default=%default]"),
  make_option(c("--top_n"), type="integer", default=10,
              help="Only show top N categories (0 for all) [default=%default]"),
  make_option(c("--title"), type="character", default="",
              help="Custom plot title [default=%default]"),
  make_option(c("--subtitle"), type="character", default="",
              help="Custom plot subtitle [default=%default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

print("Options parsed...")
print(sprintf("Input file: %s", opt$input))
print(sprintf("Output directory: %s", opt$output))

main <- function() {
  print("Entering main function...")
  
  tryCatch({
    print("Checking input file...")
    if (!file.exists(opt$input)) {
      stop(sprintf("Input file does not exist: %s", opt$input))
    }
    
    print("Reading data...")
    df <- read.csv(opt$input)
    print(sprintf("Read %d rows", nrow(df)))
    
    # Print column names and unique values
    print("Column names in dataset:")
    print(colnames(df))
    print("Unique values in model_category:")
    print(unique(df$model_category))
    print("Unique values in category (if exists):")
    if ("category" %in% colnames(df)) {
      print(unique(df$category))
    } else {
      print("No 'category' column found")
      # Try to find similar column names
      possible_category_cols <- colnames(df)[grep("category|type|class", tolower(colnames(df)))]
      print("Possible category-like columns:")
      print(possible_category_cols)
    }
    
    # Create output directory if it doesn't exist
    dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
    log_info(sprintf("Output directory: %s", opt$output))
    
    # Data preprocessing
    log_info("Preprocessing data...")
    
    # Standardize model names with more comprehensive mapping
    model_mapping <- c(
      # AlphaFold2 based
      "Alphafold2" = "Custom hallucination",
      "Alphafold2 backpropagation + alphafold2" = "Custom hallucination",
      "Alphafold2 backpropagation + alphafold2 + filter" = "Custom hallucination",
      
      # MPNN based
      "Ligand/proteinmpnn" = "ProteinMPNN",
      "Ligand/proteinmpnn/other" = "ProteinMPNN",
      "Ligand/proteinmpnn/other + alphafold2" = "ProteinMPNN",
      "Proteinmpnn + rosetta + esm2 mlde" = "ProteinMPNN",
      
      # RFdiffusion combinations
      "Rfdiffusion + ligand/proteinmpnn" = "RFdiffusion",
      "Rfdiffusion + ligand/proteinmpnn + alphafold2" = "RFdiffusion",
      "Rfdiffusion + ligand/proteinmpnn + esmfold" = "RFdiffusion",
      "Rfdiffusion + ligand/proteinmpnn + esm3" = "RFdiffusion",
      "Rfdiffusion + ligand/proteinmpnn/other + alphafold2" = "RFdiffusion",
      
      # ESM based
      "Esm2" = "ESM2",
      "Esm2 mlde" = "ESM2",
      "Esm3 conditional generation" = "ESM3",
      
      # Custom/Other methods
      "Custom active learning" = "Custom active learning",
      "Custom active learnig" = "Custom active learning",
      "Custom generative" = "Custom generative",
      "Custom diffusion" = "Custom diffusion",
      "Custom plm" = "Custom PLM",
      "Custom mlde" = "Custom MLDE",
      
      # Rosetta based
      "Rosetta" = "Rosetta",
      "Rosetta + esm2 mlde" = "Rosetta",
      
      # BindCraft
      "Bindcraft/rso" = "BindCraft/RSO",
      
      # EvoDiff
      "Evodiff conditional generation" = "EvoDiff",
      "Evodiff + esm3 + haddock + alphafold2" = "EvoDiff",
      
      # Other categories
      "Rational design" = "Rational design",
      "Other hallucination" = "Other hallucination",
      "Carbonara" = "Custom generative",
      "Raygun + protrek" = "Custom generative",
      "Lm agent" = "Other methods",
      "Timed" = "Custom generative",
      "Trained surrogate + evoprotgrad" = "Custom active learning",
      "Finetuned nach0-pc" = "Custom PLM",
      "Finetuned plm + rl" = "Custom PLM",
      "Md + docking" = "Other methods",
      
      # Additional mappings for the new entries
      "Alphafold2 backpropagation + ligand/proteinmpnn/other" = "Custom hallucination",
      "Alphafold2 backpropagation + ligand/proteinmpnn/other + alphafold2" = "Custom hallucination"
    )
    
    if ("model_category" %in% colnames(df)) {
      print("Found model_category column")
      print("Before grouping, unique values:")
      print(unique(df$model_category))
      
      df$model_names <- df$model_category  # Create new column
      for (old_name in names(model_mapping)) {
        df$model_names[tolower(df$model_names) == tolower(old_name)] <- model_mapping[old_name]
      }
      
      print("\nAfter grouping, unique values:")
      print(sort(unique(df$model_names)))
      
      # Also print counts for each group
      print("\nCounts per group:")
      print(table(df$model_names))
    }
    
    # Filter data
    df_models <- df[df$model_names != 'not mentioned' & df$model_names != 'Not mentioned',]
    
    # Only create category plot if we have the category column
    if ("category" %in% colnames(df)) {
      df_cat <- df[df$category != 'not mentioned' & df$category != 'Not mentioned',]
      print(sprintf("After filtering: %d model entries, %d category entries", 
                    nrow(df_models), nrow(df_cat)))
      
      # Create both plots
      create_and_save_plots(df_models, df_cat)
    } else {
      print(sprintf("After filtering: %d model entries (no category data available)", 
                    nrow(df_models)))
      
      # Create only model plot
      create_and_save_model_plot(df_models)
    }
    
    log_info("Script completed successfully")
    
  }, error = function(e) {
    log_error(sprintf("Error in main function: %s", e$message))
    stop(e)
  })
}

create_and_save_model_plot <- function(df_models) {
  # Define colors for rounds
  protein_colors <- c("#8B90DD", "#8CD2F4", "#3E6175", "#4A90E2")
  
  # Order by frequency
  df_models$model_names <- reorder(df_models$model_names, df_models$model_names, function(x) -length(x))
  
  # Create model plot
  bar_plot <- ggplot(df_models, aes(x = model_names, fill = factor(round))) +
    geom_bar(color = 'black') +
    geom_text(stat = 'count', aes(label = ..count..), 
              position = position_stack(vjust = 0.5), color = "white") +
    scale_fill_manual(values = protein_colors,
                     name = "Round",
                     breaks = c("1", "2", "3", "4"),
                     labels = c("Round 1", "Round 2", "Round 3", "Round 4")) +
    labs(
      title = "Model representation",
      subtitle = "Number of designs submitted by round",
      x = "Model types",
      y = "Number of designs"
    ) +
    theme_protein_ml() +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.just = "center",
      legend.margin = margin(6, 6, 6, 6),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA)
    )
  
  # Save plots
  ggsave(file.path(opt$output, 'bar_plot_submissions.png'), 
         bar_plot, 
         width = opt$width/300, 
         height = opt$height/300, 
         dpi = opt$res)
  
  ggsave(file.path(opt$output, 'bar_plot_submissions.svg'), 
         bar_plot, 
         width = opt$width/300, 
         height = opt$height/300)
}

create_bar_plot <- function(data, x_col, title) {
  ggplot(data, aes_string(x = x_col, fill = "selected")) +
    geom_bar(color = 'black') +
    geom_text(stat = 'count', aes(label = ..count..), 
              position = position_stack(vjust = 0.5), color = "white") +
    scale_fill_manual(values = protein_colors,
                     name = "Selection status",
                     breaks = c("Top 100", "Adaptyv selection", "No"),
                     labels = c("Top 100", "Adaptyv selection", "Not selected")) +
    scale_y_continuous() +
    labs(
      title = title,
      subtitle = "Number of designs submitted vs. selected for validation",
      x = ifelse(x_col == "model_names", "Model types", "Design category"),
      y = "Number of designs"
    ) +
    theme_protein_ml() +
    theme(
      legend.position = "right",
      legend.justification = "center",
      legend.box.just = "center",
      legend.margin = margin(6, 6, 6, 6),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA)
    ) +
    guides(fill = guide_legend(title = "Selection status"))
}

save_plots <- function(bar_plot, bar_plot_cat, output_dir, height, width, res) {
  # Save PNG files
  png(file.path(output_dir, 'bar_plot_submissions.png'), 
      bg = 'white', height = height, width = width, res = res)
  print(bar_plot)
  dev.off()
  
  png(file.path(output_dir, 'bar_plot_cat.png'), 
      bg = 'white', height = height, width = width, res = res)
  print(bar_plot_cat)
  dev.off()
  
  # Save SVG files
  svg(file.path(output_dir, 'bar_plot_submissions.svg'), 
      bg = 'white', height = height, width = width)
  print(bar_plot)
  dev.off()
  
  svg(file.path(output_dir, 'bar_plot_cat.svg'), 
      bg = 'white', height = height, width = width)
  print(bar_plot_cat)
  dev.off()
}

# Add print statement before running main
print("About to run main function...")
if (sys.nframe() == 0) {
  main()
}
print("Script completed.")