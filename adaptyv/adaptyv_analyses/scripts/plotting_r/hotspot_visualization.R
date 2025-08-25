#!/usr/bin/env Rscript

library(data.table)
library(reticulate)
library(optparse)
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(stringr)
library(sysfonts)
library(showtext)

# Source the Adaptyv theme
source("scripts/plotting/blog_post_theme.R")

# Initialize Python environment from conda
use_condaenv("/opt/homebrew/Caskroom/miniforge/base/envs/scicomp")

# Import required modules
matplotlib <- import("matplotlib")
plt <- import("matplotlib.pyplot")
np <- import("numpy")
scipy_special <- import("scipy.special")

# Add Roboto font
font_add_google("Roboto", "Roboto")
showtext_auto()

# Define domain boundaries and colors with stronger base coloring
domains <- list(
    I = list(range=1:165, color="#3E6175"),    # Dark teal
    II = list(range=166:310, color="#56A6D4"),  # Light blue
    III = list(range=311:480, color="#8B90DD"), # Periwinkle
    IV = list(range=481:620, color="#2A4DD0")   # Dark blue
)

option_list <- list(
    make_option(c("--input"), type="character", 
                default="/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/processed/all_submissions.csv",
                help="Input CSV file path"),
    make_option(c("--structure"), type="character", 
                default="./data/raw/structures/002_2024/02f25cbe-186c-4a52-8985-05e13f4d183b_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_000.pdb",
                help="Input PDB structure file"),
    make_option(c("--round"), type="character", default=2,
                help="Filter by round (1, 2, or both)"),
    make_option(c("--binding"), type="character", default="all",
                help="Filter by binding (Yes, No, or all)"),
    make_option(c("--binding_strength"), type="character", default="all",
                help="Filter by binding strength (Weak, Medium, Strong, or all)"),
    make_option(c("--title"), type="character", default="",
                help="Custom title for the plot"),
    make_option(c("--subtitle"), type="character", default="",
                help="Custom subtitle for the plot"),
    make_option(c("--format"), type="character", default="png",
                help="Output format (png or svg)"),
    make_option(c("--width"), type="integer", default=3600,
                help="Plot width in pixels"),
    make_option(c("--height"), type="integer", default=3600,
                help="Plot height in pixels"),
    make_option(c("--res"), type="integer", default=600,
                help="Plot resolution in DPI")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

main <- function() {
    # Read data
    data <- fread(opt$input)
    
    # Function to process data for a specific filter condition
    process_data <- function(data, binding_filter=NULL, strength_filter=NULL, de_novo_filter=NULL, group_name) {
        filtered_data <- data
        if (!is.null(binding_filter)) {
            filtered_data <- filtered_data[binding %in% binding_filter]
        }
        if (!is.null(strength_filter)) {
            filtered_data <- filtered_data[binding_strength %in% strength_filter]
        }
        if (!is.null(de_novo_filter)) {
            filtered_data <- filtered_data[de_novo %in% de_novo_filter]
        }
        
        domain_counts <- sapply(names(domains), function(x) 0)
        names(domain_counts) <- names(domains)
        
        for(entry in filtered_data$target_residue_list) {
            if(is.character(entry)) {
                nums <- as.numeric(trimws(strsplit(gsub("\\[|\\]", "", entry), ",")[[1]]))
            } else {
                nums <- entry
            }
            
            for(domain_name in names(domains)) {
                domain <- domains[[domain_name]]
                count_in_domain <- sum(nums %in% domain$range)
                domain_counts[domain_name] <- domain_counts[domain_name] + count_in_domain
            }
        }
        
        total_hotspots <- sum(domain_counts)
        
        plot_df <- data.frame(
            domain = sprintf("Domain %s", names(domain_counts)),
            count = as.numeric(domain_counts),
            percentage = round((as.numeric(domain_counts) / total_hotspots) * 100, 1),
            group = group_name,
            color = sapply(names(domain_counts), function(x) domains[[x]]$color)
        )
        return(plot_df)
    }
    
    # Process data for each condition
    strong_binders <- process_data(data, binding_filter="Yes", 
                                 strength_filter="Strong", group_name="Strong binders")
    all_binders <- process_data(data, binding_filter="Yes", 
                              group_name="All binders")
    all_tested <- process_data(data, binding_filter=c("Yes", "No", "Unknown"), 
                             group_name="All tested")
    all_submissions <- process_data(data, group_name="All submissions")
    de_novo_binders <- process_data(data, binding_filter="Yes", 
                                  de_novo_filter="De novo", group_name="De novo binders")
    
    # Combine all data
    plot_df <- rbind(strong_binders, all_binders, all_tested, all_submissions, de_novo_binders)
    
    # Set factor levels for ordering
    plot_df$group <- factor(plot_df$group, 
                          levels=c("All submissions", "All tested", "All binders", "Strong binders", "De novo binders"))
    
    # Create stacked barplot
    p <- ggplot(plot_df, aes(x=group, y=percentage, fill=domain)) +
        geom_bar(stat="identity", position="stack", color="black", linewidth=0.5) +
        geom_text(data = subset(plot_df, percentage > 0),
                 aes(label=sprintf("%.1f%%", percentage)), 
                 position=position_stack(vjust=0.5),
                 color="white", size=3.5) +
        scale_fill_manual(values=unique(plot_df$color)) +
        labs(title = "Domain distribution of targeted sites",
             subtitle = "Indicating the EGFR domain bias for all binders",
             x="", y="Percentage of hotspots", fill="Domain") +
        adaptyv_theme() +
        theme(
            axis.text.x = element_text(angle=45, hjust=1),
            legend.position="right"
        )
    
    # Save plot
    dir.create("plots/hotspots", showWarnings = FALSE, recursive = TRUE)
    
    filename <- file.path("plots/hotspots", 
                         sprintf("domain_distribution_stacked.%s", opt$format))
    
    if (opt$format == "png") {
        png(filename, bg = 'white', 
            height = opt$height, width = opt$width, res = opt$res)
        print(p)
        dev.off()
        message(sprintf('Saved PNG plot to %s', filename))
    } else if (opt$format == "svg") {
        svg(filename, bg = 'transparent')
        print(p)
        dev.off()
        message(sprintf('Saved SVG plot to %s', filename))
    }
}

if (sys.nframe() == 0) {
    main()
}