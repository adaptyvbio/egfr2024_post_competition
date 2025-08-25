#!/usr/bin/env Rscript

library(data.table)
library(reticulate)
library(optparse)

source("./scripts/plotting/blog_post_theme.R")

use_condaenv("/opt/homebrew/Caskroom/miniforge/base/envs/scicomp")

matplotlib <- import("matplotlib")
plt <- import("matplotlib.pyplot")
np <- import("numpy")
scipy_special <- import("scipy.special")

domains <- list(
    I = list(range=1:165, color="#3E6175"),
    II = list(range=166:310, color="#56A6D4"),
    III = list(range=311:480, color="#8B90DD"),
    IV = list(range=481:620, color="#2A4DD0")
)

py_run_string('
def plot_pseudo_3D(xyz, colors, ax, line_w=2.0):
    import numpy as np
    import matplotlib.patheffects as patheffects
    from matplotlib import collections as mcoll
    from matplotlib.colors import to_rgba_array
    
    xyz = np.asarray(xyz)
    colors = to_rgba_array(colors)
    
    segments = np.zeros((len(xyz)-1, 2, 3))
    segments[:,0] = xyz[:-1]
    segments[:,1] = xyz[1:]
    
    segments_xy = segments[...,:2]
    segments_z = segments[...,2].mean(axis=1)
    
    order = segments_z.argsort()
    
    segment_colors = np.zeros((len(xyz)-1, 4))
    segment_colors[...,:3] = (colors[:-1,:3] + colors[1:,:3]) / 2
    segment_colors[...,3] = 1.0
    
    z_norm = (segments_z - segments_z.min()) / (segments_z.max() - segments_z.min())
    segment_colors[...,:3] = segment_colors[...,:3] * (0.7 + 0.3 * z_norm[:,None])
    
    lines = mcoll.LineCollection(
        segments_xy[order],
        colors=segment_colors[order],
        linewidths=line_w,
        path_effects=[
            patheffects.Stroke(linewidth=line_w+1, foreground="white"),
            patheffects.Normal()
        ]
    )
    
    ax.add_collection(lines)
    
    margin = line_w * 2
    ax.set_xlim(xyz[:,0].min() - margin, xyz[:,0].max() + margin)
    ax.set_ylim(xyz[:,1].min() - margin, xyz[:,1].max() + margin)
    
    return lines
')

option_list <- list(
    make_option(c("--input"), type="character", 
                default="/Users/tudorcotet/Documents/Adaptyv/competition_paper/data/processed/all_submissions.csv",
                help="Input CSV file path"),
    make_option(c("--structure"), type="character", 
                default="./data/raw/structures/002_2024/02f25cbe-186c-4a52-8985-05e13f4d183b_unrelaxed_rank_001_alphafold2_multimer_v3_model_3_seed_000.pdb",
                help="Input PDB structure file"),
    make_option(c("--round"), type="character", default=2,
                help="Filter by round (1, 2, or both)"),
    make_option(c("--binding"), type="character", default="Yes",
                help="Filter by binding (Yes, No, or all)"),
    make_option(c("--binding_strength"), type="character", default="all",
                help="Filter by binding strength (Weak, Medium, Strong, or all)"),
    make_option(c("--de_novo"), type="character", default="Existing binder",
                help="Filter by de novo status (De novo, Existing binder, or all)"),
    make_option(c("--title"), type="character", default="Targeted site distribution on EGFR",
                help="Custom title for the plot"),
    make_option(c("--subtitle"), type="character", default="For existing binders",
                help="Custom subtitle for the plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

interpolate_color <- function(color1, color2, weight) {
    rgb1 <- col2rgb(color1)
    rgb2 <- col2rgb(color2)
    
    rgb_interp <- round(rgb1 * weight + rgb2 * (1 - weight))
    
    rgb_interp <- rgb(rgb_interp[1], rgb_interp[2], rgb_interp[3])
    return(rgb_interp)
}

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

main <- function() {
    data <- fread(opt$input)
    message("Initial data dimensions: ", paste(dim(data), collapse=" x "))
    
    if (opt$round != "both") {
        round_num <- as.numeric(opt$round)
        data <- data[round == round_num]
        message(sprintf("Filtered for round %d", round_num))
    }
    
    if (opt$binding != "all") {
        data <- data[binding == opt$binding]
        message(sprintf("Filtered for binding = %s", opt$binding))
    }
    
    if (opt$binding_strength != "all") {
        data <- data[binding_strength == opt$binding_strength]
        message(sprintf("Filtered for binding_strength = %s", opt$binding_strength))
    }
    
    if (opt$de_novo != "all") {
        data <- data[de_novo == opt$de_novo]
        message(sprintf("Filtered for de_novo = %s", opt$de_novo))
    }
    
    message("Final data dimensions: ", paste(dim(data), collapse=" x "))
    
    if(nrow(data) > 0 && "target_residue_list" %in% names(data)) {
        residue_lists <- lapply(data$target_residue_list, function(x) {
            if(is.character(x)) {
                nums <- strsplit(gsub("\\[|\\]", "", x), ",")[[1]]
                as.numeric(trimws(nums))
            } else {
                x
            }
        })
        
        all_residues <- unlist(residue_lists)
        residue_freq <- table(all_residues)
        
        hotspot_freq <- as.numeric(residue_freq) / max(as.numeric(residue_freq))
        names(hotspot_freq) <- names(residue_freq)
        
        message(sprintf("Found %d unique residues", length(unique(all_residues))))
        message(sprintf("Maximum frequency: %d occurrences", max(as.numeric(residue_freq))))
    } else {
        message("No target_residue_list column found")
        hotspot_freq <- numeric(0)
    }
    
    pdb_lines <- readLines(opt$structure)
    xyz <- np$array(matrix(0, nrow=length(pdb_lines), ncol=3))
    
    ca_idx <- 1
    residue_numbers <- c()
    for(line in pdb_lines) {
        if(startsWith(line, "ATOM") && grepl("CA", line)) {
            xyz[ca_idx,] <- as.numeric(c(
                substr(line, 31, 38),
                substr(line, 39, 46),
                substr(line, 47, 54)
            ))
            residue_numbers[ca_idx] <- as.numeric(substr(line, 23, 26))
            ca_idx <- ca_idx + 1
        }
    }
    xyz <- xyz[1:(ca_idx-1),]
    residue_numbers <- residue_numbers[1:(ca_idx-1)]
    
    base_gray <- "#E8E8E8"
    colors <- rep(base_gray, nrow(xyz))
    
    for(res_num in residue_numbers) {
        for(domain_name in names(domains)) {
            domain <- domains[[domain_name]]
            if(res_num %in% domain$range) {
                if(as.character(res_num) %in% names(hotspot_freq)) {
                    freq <- hotspot_freq[as.character(res_num)]
                    
                    rgb1 <- col2rgb(base_gray)
                    rgb2 <- col2rgb(domain$color)
                    
                    mixed_rgb <- rgb1 * (1 - freq) + rgb2 * freq
                    mixed_color <- rgb(mixed_rgb[1]/255, mixed_rgb[2]/255, mixed_rgb[3]/255)
                    
                    colors[residue_numbers == res_num] <- mixed_color
                }
                break
            }
        }
    }
    
    plt$figure(figsize=c(10, 9))
    ax <- plt$gca()
    ax$set_aspect('equal')
    ax$set_facecolor("white")
    
    if (opt$title == "") {
        title <- "Hotspot analysis"
    } else {
        title <- opt$title
    }

    if (opt$subtitle == "") {
        subtitle <- sprintf("Round: %s | Binding: %s | Strength: %s | De novo: %s", 
                          opt$round, opt$binding, opt$binding_strength, opt$de_novo)
    } else {
        subtitle <- opt$subtitle
    }
    
    plt$suptitle(title, y=0.95, fontsize=16, fontweight="bold")
    plt$title(subtitle, color="#3D7A9A", fontsize=12, pad=20)
    
    ax$set_xticks(np$array(numeric()))
    ax$set_yticks(np$array(numeric()))
    
    tryCatch({
        py$plot_pseudo_3D(
            xyz = np$array(xyz),
            colors = colors,
            ax = ax,
            line_w = 2.5
        )
        
        for(domain_name in names(domains)) {
            domain <- domains[[domain_name]]
            domain_range <- domain$range
            domain_xyz <- xyz[residue_numbers %in% domain_range,]
            if(nrow(domain_xyz) > 0) {
                center <- colMeans(domain_xyz)
                ax$text(
                    center[1], center[2],
                    sprintf("Domain %s", domain_name),
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=12,
                    fontweight="bold",
                    color=domain$color,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none")
                )
            }
        }
        
        dir.create("plots/hotspots", showWarnings = FALSE, recursive = TRUE)
        plt$savefig("plots/hotspots/structure.png", dpi=300, bbox_inches="tight", 
                   facecolor="white", edgecolor="none")
        message("Saved visualization to plots/hotspots/structure.png")
    }, error = function(e) {
        message("Error in plotting: ", e$message)
        print(py_last_error())
    })
}

if (sys.nframe() == 0) {
    main()
}