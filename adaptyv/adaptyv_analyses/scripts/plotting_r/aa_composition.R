#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(jsonlite)

# Import theme and colors
font_add_google("Roboto", "Roboto")
showtext_auto()

# Define color palette for de novo vs existing binders vs non-binders (matching violin_properties.R)
binding_status_colors <- c(
  "De novo" = "#67AFD8",
  "Existing binder" = "#8D92DE",
  "Non-binder" = "#E5E7EB"
)

# Define amino acid classes and their colors
aa_classes <- list(
  "Hydrophobic" = c("A", "I", "L", "M", "V", "F", "W", "Y"),  # Nonpolar (aliphatic + aromatic)
  "Polar" = c("S", "T", "N", "Q"),  # Polar uncharged
  "Basic" = c("K", "R", "H"),  # Positively charged at pH 7.4
  "Acidic" = c("D", "E"),  # Negatively charged at pH 7.4
  "Special" = c("G", "P", "C")  # Glycine (no side chain), Proline (cyclic), Cysteine (disulfide bonds)
)

aa_class_colors <- c(
  "Hydrophobic" = "#2A4DD0",
  "Polar" = "#A7C1FB",
  "Basic" = "#C4B2FB",
  "Acidic" = "#DC7A73",
  "Special" = "#69B7A7"
)

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_markdown(size = 16, color = "#3D7A9A", hjust = 0.5, margin = margin(b = 20)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 0.5),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      legend.position = "right",
      legend.box.background = element_rect(color = "black"),
      plot.margin = margin(20, 20, 20, 20),
      strip.text = element_text(size = 16, face = "bold", color = "#3D7A9A"),
      strip.background = element_rect(fill = "white", color = NA),
      panel.spacing = unit(2, "lines")
    )
}

# Command line options
option_list <- list(
  make_option(c("--input"), type="character", default="./data/processed/all_submissions.csv",
              help="Input CSV file path"),
  make_option(c("--output"), type="character", default="./plots/aa_composition/",
              help="Output directory path"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=8400,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=3200,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=600,
              help="Plot resolution (DPI)"),
  make_option(c("--round"), type="character", default='both',
              help="Filter by round (1, 2, or both)"),
  make_option(c("--title"), type="character", default="Amino acid composition comparison",
              help="Custom plot title"),
  make_option(c("--subtitle"), type="character", default="Distribution across de novo and existing binders",
              help="Custom plot subtitle")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to parse binder residues from JSON-like string
parse_binder_residues <- function(json_str) {
  if (is.na(json_str) || json_str == "") return(NULL)
  
  tryCatch({
    # Extract amino acids using regex - looking for single quotes between colons and commas/braces
    matches <- gregexpr(":\\s*'([A-Z])'", json_str)
    amino_acids <- regmatches(json_str, matches)[[1]]
    # Clean up to get just the amino acids
    amino_acids <- gsub(":\\s*'|'", "", amino_acids)
    return(amino_acids)
  }, error = function(e) {
    log_warn(paste("Failed to parse binder residues:", json_str))
    return(NULL)
  })
}

# Function to count amino acids in a sequence
count_sequence_aas <- function(sequence) {
  if (is.na(sequence) || sequence == "") return(NULL)
  
  # Split sequence into individual amino acids
  amino_acids <- strsplit(sequence, "")[[1]]
  return(amino_acids)
}

# Function to get amino acid class
get_aa_class <- function(aa) {
  for (class_name in names(aa_classes)) {
    if (aa %in% aa_classes[[class_name]]) {
      return(class_name)
    }
  }
  return("Other")
}

main <- function() {
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
  
  # Parse both full sequence and interface residues and count amino acids
  full_seq_counts <- data.table()
  interface_counts <- data.table()
  
  for (i in 1:nrow(df)) {
    # Count full sequence amino acids
    full_seq_aas <- count_sequence_aas(df$sequence[i])
    if (!is.null(full_seq_aas)) {
      counts <- table(full_seq_aas)
      aa_dt <- data.table(
        amino_acid = names(counts),
        count = as.numeric(counts),
        design_category = df$design_category[i],
        type = "Full sequence"
      )
      full_seq_counts <- rbind(full_seq_counts, aa_dt)
    }
    
    # Count interface residues
    interface_aas <- parse_binder_residues(df$binder_residues[i])
    if (!is.null(interface_aas)) {
      counts <- table(interface_aas)
      aa_dt <- data.table(
        amino_acid = names(counts),
        count = as.numeric(counts),
        design_category = df$design_category[i],
        type = "Interface"
      )
      interface_counts <- rbind(interface_counts, aa_dt)
    }
  }
  
  # Combine both counts
  aa_counts <- rbind(full_seq_counts, interface_counts)
  
  # Add amino acid class information
  aa_counts[, aa_class := sapply(amino_acid, get_aa_class)]
  
  # Aggregate counts by amino acid, design category, and type
  aa_summary <- aa_counts[, .(total_count = sum(count)), 
                         by = .(amino_acid, design_category, type)]
  
  # Calculate percentages within each group and type
  aa_summary[, percentage := total_count / sum(total_count) * 100, 
             by = .(design_category, type)]
  
  # Calculate standard deviation of percentages for each amino acid and type
  aa_summary[, percentage_sd := sd(percentage), by = .(amino_acid, type)]
  
  # Calculate differences between categories for each amino acid and type
  aa_summary[, max_diff := max(percentage) - min(percentage), by = .(amino_acid, type)]
  
  # Calculate pairwise differences between categories
  aa_summary[, binder_vs_nonbinder_diff := abs(max(percentage[design_category %in% c("De novo", "Existing binder")]) - 
                                              percentage[design_category == "Non-binder"]), 
            by = .(amino_acid, type)]
  
  aa_summary[, denovo_vs_existing_diff := abs(percentage[design_category == "De novo"] - 
                                             percentage[design_category == "Existing binder"]), 
            by = .(amino_acid, type)]
  
  # Calculate mean percentage for each amino acid and type
  aa_summary[, mean_pct := mean(percentage), by = .(amino_acid, type)]
  
  # Flag amino acids with significant differences (relative to their mean)
  # Now considering both overall differences and pairwise differences
  aa_summary[, is_significant := (max_diff >= 0.15 * mean_pct) |  # 15% overall difference
                                (binder_vs_nonbinder_diff >= 0.15 * mean_pct) |  # 15% binder vs non-binder
                                (denovo_vs_existing_diff >= 0.15 * mean_pct)]  # 15% de novo vs existing
  
  # Set factor levels to control stacking order (non-binders at bottom)
  aa_summary[, design_category := factor(design_category, 
                                       levels = c("Non-binder", "Existing binder", "De novo"))]
  
  # Aggregate by amino acid class
  class_summary <- aa_counts[, .(total_count = sum(count)), 
                           by = .(aa_class, design_category, type)]
  
  # Calculate percentages for classes
  class_summary[, percentage := total_count / sum(total_count) * 100, 
               by = .(design_category, type)]
  
  # Calculate standard deviation of percentages for each class
  class_summary[, percentage_sd := sd(percentage), by = .(aa_class, type)]
  
  # Calculate differences between categories for each class and type
  class_summary[, max_diff := max(percentage) - min(percentage), by = .(aa_class, type)]
  
  # Calculate mean percentage for each class and type
  class_summary[, mean_pct := mean(percentage), by = .(aa_class, type)]
  
  # Flag classes with significant differences (relative to their mean)
  class_summary[, is_significant := max_diff >= 0.15 * mean_pct]  # 15% relative difference threshold
  
  # Set factor levels for classes
  class_summary[, aa_class := factor(aa_class, levels = names(aa_classes))]
  class_summary[, design_category := factor(design_category, 
                                          levels = c("Non-binder", "Existing binder", "De novo"))]
  
  # Print analysis of percentages
  cat("\nAnalysis of amino acid compositions:\n")
  cat("=====================================\n\n")
  
  # Function to format percentage differences
  format_diff <- function(diff) {
    if(diff > 0) return(sprintf("+%.1f%%", diff))
    return(sprintf("%.1f%%", diff))
  }
  
  # Analyze each type separately
  for(current_type in unique(aa_summary$type)) {
    cat(sprintf("\n%s:\n", current_type))
    cat("----------------------------------------\n")
    
    # Get data for current type
    type_data <- aa_summary[type == current_type]
    
    # For each amino acid with significant differences
    for(aa in unique(type_data[is_significant == TRUE]$amino_acid)) {
      aa_data <- type_data[amino_acid == aa]
      
      # Get percentages for each category
      denovo_pct <- aa_data[design_category == "De novo"]$percentage
      existing_pct <- aa_data[design_category == "Existing binder"]$percentage
      nonbinder_pct <- aa_data[design_category == "Non-binder"]$percentage
      
      # Calculate differences
      binder_vs_non <- mean(c(denovo_pct, existing_pct)) - nonbinder_pct
      denovo_vs_existing <- denovo_pct - existing_pct
      
      cat(sprintf("\nAmino acid %s:\n", aa))
      cat(sprintf("  De novo: %.1f%%\n", denovo_pct))
      cat(sprintf("  Existing: %.1f%%\n", existing_pct))
      cat(sprintf("  Non-binder: %.1f%%\n", nonbinder_pct))
      cat(sprintf("  Differences:\n"))
      cat(sprintf("    Binders vs Non-binders: %s\n", format_diff(binder_vs_non)))
      cat(sprintf("    De novo vs Existing: %s\n", format_diff(denovo_vs_existing)))
    }
  }
  
  cat("\nAnalysis of amino acid type compositions:\n")
  cat("=========================================\n\n")
  
  # Analyze each type separately for classes
  for(current_type in unique(class_summary$type)) {
    cat(sprintf("\n%s:\n", current_type))
    cat("----------------------------------------\n")
    
    # Get data for current type
    type_data <- class_summary[type == current_type]
    
    # For each class
    for(cls in unique(type_data$aa_class)) {
      class_data <- type_data[aa_class == cls]
      
      # Get percentages for each category
      denovo_pct <- class_data[design_category == "De novo"]$percentage
      existing_pct <- class_data[design_category == "Existing binder"]$percentage
      nonbinder_pct <- class_data[design_category == "Non-binder"]$percentage
      
      # Calculate differences
      binder_vs_non <- mean(c(denovo_pct, existing_pct)) - nonbinder_pct
      denovo_vs_existing <- denovo_pct - existing_pct
      
      cat(sprintf("\n%s:\n", cls))
      cat(sprintf("  De novo: %.1f%%\n", denovo_pct))
      cat(sprintf("  Existing: %.1f%%\n", existing_pct))
      cat(sprintf("  Non-binder: %.1f%%\n", nonbinder_pct))
      cat(sprintf("  Differences:\n"))
      cat(sprintf("    Binders vs Non-binders: %s\n", format_diff(binder_vs_non)))
      cat(sprintf("    De novo vs Existing: %s\n", format_diff(denovo_vs_existing)))
    }
  }

  # Create amino acid composition plot
  p1 <- ggplot(aa_summary, aes(x = amino_acid, y = percentage, fill = design_category)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    geom_text(data = aa_summary[is_significant == TRUE],
              aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5),
              size = 3.5,
              color = ifelse(aa_summary[is_significant == TRUE]$design_category == "Non-binder", "black", "white"),
              fontface = "bold") +
    scale_fill_manual(values = binding_status_colors,
                     name = "Binder type") +
    facet_wrap(~type, ncol = 2, scales = "free_y",
               labeller = as_labeller(c(
                 "Full sequence" = "Full sequence composition",
                 "Interface" = "Interface composition"
               ))) +
    labs(title = "Amino acid composition",
         subtitle = "Distribution comparison between full sequence and interface residues",
         x = "Amino acid",
         y = "Percentage (%)") +
    adaptyv_theme() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Save amino acid composition plot
  filename1 <- file.path(opt$output, 
                       sprintf("aa_composition_round%s.%s", 
                             opt$round,
                             opt$format))
  
  if (opt$format == "png") {
    png(filename1, width = opt$width, height = opt$height, res = opt$res)
    print(p1)
    dev.off()
  } else {
    svg(filename1, width = opt$width/opt$res, height = opt$height/opt$res)
    print(p1)
    dev.off()
  }
  
  # Create amino acid class plot
  p2 <- ggplot(class_summary, aes(x = aa_class, y = percentage, fill = design_category)) +
    geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
    geom_text(data = class_summary[percentage >= 2],
              aes(label = sprintf("%.1f%%", percentage)),
              position = position_stack(vjust = 0.5),
              size = 3,
              color = ifelse(class_summary[percentage >= 2]$design_category == "Non-binder", "black", "white"),
              fontface = "bold") +
    scale_fill_manual(values = binding_status_colors,
                     name = "Binder type") +
    facet_wrap(~type, ncol = 2, scales = "free_y",
               labeller = as_labeller(c(
                 "Full sequence" = "Full sequence composition",
                 "Interface" = "Interface composition"
               ))) +
    labs(title = "Amino acid type composition",
         subtitle = "Distribution of amino acid types in full sequence and interface residues",
         x = "Amino acid type",
         y = "Percentage (%)") +
    adaptyv_theme() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Save amino acid class plot
  filename2 <- file.path(opt$output, 
                       sprintf("aa_type_composition_round%s.%s", 
                             opt$round,
                             opt$format))
  
  if (opt$format == "png") {
    png(filename2, width = opt$width, height = opt$height, res = opt$res)
    print(p2)
    dev.off()
  } else {
    svg(filename2, width = opt$width/opt$res, height = opt$height/opt$res)
    print(p2)
    dev.off()
  }
  
  # Calculate enrichment ratios
  enrichment_data <- merge(
    class_summary[type == "Interface", .(aa_class, design_category, interface_pct = percentage)],
    class_summary[type == "Full sequence", .(aa_class, design_category, full_seq_pct = percentage)],
    by = c("aa_class", "design_category")
  )
  
  enrichment_data[, enrichment_ratio := interface_pct / full_seq_pct]
  
  # Make enrichment plot more compact but ensure text visibility
  opt$width <- 5000  # Increased from 4000
  opt$height <- 3000  # Increased from 2400
  
  # Create enrichment ratio plot with adjusted design
  p3 <- ggplot(enrichment_data, aes(x = aa_class, y = enrichment_ratio, fill = design_category)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = sprintf("%.1fx", enrichment_ratio)),
              position = position_dodge(width = 0.9),
              vjust = -0.2,
              size = 2.5,
              fontface = "bold") +
    scale_fill_manual(values = binding_status_colors,
                     name = "Binder type") +
    labs(title = "Interface enrichment",
         subtitle = "Ratio of amino acid type frequencies at interface vs full sequence",
         x = "Amino acid type",
         y = "Enrichment ratio (interface / full sequence)") +
    adaptyv_theme() +
    theme(
      plot.margin = margin(30, 20, 20, 20),  # Increased top margin for title
      panel.spacing = unit(1.5, "lines"),  # Slightly increased panel spacing
      axis.title.y = element_text(margin = margin(r = 10))  # Add margin to y-axis title
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))  # Increased top padding
  
  # Save enrichment plot
  filename3 <- file.path(opt$output, 
                        sprintf("aa_type_enrichment_round%s.%s", 
                              opt$round,
                              opt$format))
  
  if (opt$format == "png") {
    png(filename3, width = opt$width, height = opt$height, res = opt$res)
    print(p3)
    dev.off()
  } else {
    svg(filename3, width = opt$width/opt$res, height = opt$height/opt$res)
    print(p3)
    dev.off()
  }
  
  log_info(sprintf("Plots saved to %s, %s, and %s", filename1, filename2, filename3))
}

if (sys.nframe() == 0) {
  main()
}
