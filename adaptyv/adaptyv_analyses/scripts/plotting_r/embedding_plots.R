#!/usr/bin/env Rscript

library(ggplot2)
library(data.table)
library(showtext)
library(scales)
library(ggtext)
library(dplyr)
library(logger)
library(optparse)
library(umap)
library(Rtsne)
library(reticulate)
library(rhdf5)
library(openssl)  # for sha256 hashing
library(ggforce)  # for convex hulls and advanced plotting
library(future)   # for parallel processing
library(furrr)    # for parallel map operations
library(tictoc)   # for timing
library(arrow)    # for parquet files
library(stringdist)  # for sequence similarity
use_condaenv("scicomp", required = TRUE)
pacmap <- import("pacmap")

# Set up parallel processing with memory limits
plan(multisession, workers = min(4, parallel::detectCores() - 1))  # Use fewer workers to avoid memory issues

# Set up logging
log_threshold(DEBUG)

font_add_google("Roboto", "Roboto")
showtext_auto()

adaptyv_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Roboto", color = "#333333"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
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
      legend.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      #legend.box.background = element_rect(color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )
}

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
  "None" = '#E5E7EB',
  "Low" = "#9EA2AF",
  "Medium" = "#E9C435",
  "High" = "#69B7A7"
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
  "Hallucination" = "#3E6175"
)

de_novo_colors <- c(
  "De novo" = "#56A6D4",
  "Existing binder" = "#FFB347"
)

metric_colors <- c("ESM2 PLL" = "#8CD2F4", "ipTM" = "#8B90DD", "iPAE" = "#3E6175")

adaptyv_colors <- unique(unlist(c(binding_colors, expression_colors, selection_status_colors, design_category_colors, metric_colors)))

homology_colors <- c(
  "EGF-like" = "#FFB6C1",      # Light pink
  "TGFα-like" = "#B0E0E6",     # Light blue
  "Antibody-like" = "#FFE4B5"  # Light orange/peach
)

#keep like this
round_colors <- c(
  "2" = "#8CD2F4",  # blue
  "1" = "#8B90DD"   # purple
)

reference_sequences <- list(
  egf = "NSDSECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCYRDLKWWELR",
  tgf = "VVSHFNDCPDSHTQFCFHGTCRFLVQEDKPACVCHSGYVGARCEHADLLA",
  cetuximab = paste0(
    # VH only
    "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYN",
    "TPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAA"
  )
)

calculate_sequence_similarity <- function(sequence, reference) {
  # Convert sequences to uppercase for consistency
  sequence <- toupper(sequence)
  reference <- toupper(reference)
  
  # Initialize the scoring matrix
  n <- nchar(sequence)
  m <- nchar(reference)
  score_matrix <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  # Gap penalty and match/mismatch scores
  gap_penalty <- -2
  match_score <- 1
  mismatch_score <- -1
  
  # Split sequences into character vectors
  seq_chars <- strsplit(sequence, "")[[1]]
  ref_chars <- strsplit(reference, "")[[1]]
  
  # Fill the first row and column (gap penalties)
  score_matrix[1,] <- 0:(m) * gap_penalty
  score_matrix[,1] <- 0:(n) * gap_penalty
  
  # Fill the scoring matrix
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      match <- score_matrix[i-1, j-1] + 
        ifelse(seq_chars[i-1] == ref_chars[j-1], match_score, mismatch_score)
      delete <- score_matrix[i-1, j] + gap_penalty
      insert <- score_matrix[i, j-1] + gap_penalty
      score_matrix[i, j] <- max(match, delete, insert)
    }
  }
  
  # Calculate maximum possible score for normalization
  max_possible_score <- min(n, m) * match_score
  
  # Return normalized similarity score
  similarity <- score_matrix[n + 1, m + 1] / max_possible_score
  return(max(0, similarity))  # Ensure non-negative similarity
}

# Command line options
option_list <- list(
  make_option(c("--embeddings"), type="character", default='./data/processed/embeddings/saprot.parquet',
              help="Input parquet/csv file with sequence and embedding columns"),
  make_option(c("--metadata"), type="character", default="./data/processed/all_submissions.csv",
              help="Metadata CSV file path"),
  make_option(c("--round"), type="character", default='both',
              help="Filter by round (1, 2, or both)"),
  make_option(c("--output"), type="character", default="./plots/embeddings/",
              help="Output directory path"),
  make_option(c("--color_by"), type="character", default='round',
              help="Column to color points by"),
  make_option(c("--method"), type="character", default="umap",
              help="Dimensionality reduction method (pca, umap, tsne, or pacmap)"),
  make_option(c("--title"), type="character", default=NULL,
              help="Custom plot title"),
  make_option(c("--subtitle"), type="character", default=NULL,
              help="Custom plot subtitle"),
  make_option(c("--format"), type="character", default="svg",
              help="Output format (png or svg)"),
  make_option(c("--width"), type="integer", default=1600,
              help="Plot width in pixels"),
  make_option(c("--height"), type="integer", default=1200,
              help="Plot height in pixels"),
  make_option(c("--res"), type="integer", default=300,
              help="Plot resolution (DPI)"),
  make_option(c("--tested_only"), type="logical", default=FALSE,
              help="Filter for tested designs only (Adaptyv selection or Top 100) [default=%default]"),
  make_option(c("--embedding_type"), type="character", default="saprot",
              help="Type of embedding (esm2, esm1b, etc)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to perform dimensionality reduction with memory management
reduce_dimensions <- function(embeddings, method) {
  log_info(sprintf("Starting %s reduction", method))
  
  # Convert embeddings to numeric matrix and check for validity
  embeddings_matrix <- tryCatch({
    # If embeddings is already a matrix, use it directly
    if (is.matrix(embeddings)) {
      log_info(sprintf("Using provided matrix of size %d x %d", nrow(embeddings), ncol(embeddings)))
      embeddings
    } else {
      # Otherwise, convert to matrix
      matrix <- as.matrix(embeddings)
      log_info(sprintf("Created matrix of size %d x %d", nrow(matrix), ncol(matrix)))
      matrix
    }
  }, error = function(e) {
    log_error("Failed to convert embeddings to numeric matrix")
    log_error(e$message)
    stop("Embeddings must contain only numeric values")
  })
  
  # Check for NAs or non-finite values
  if (any(is.na(embeddings_matrix)) || any(!is.finite(embeddings_matrix))) {
    log_error("Embeddings contain NA or non-finite values")
    stop("Invalid values found in embeddings")
  }

  result <- tryCatch({
  if (method == "pca") {
      log_info("Computing PCA")
    pca <- prcomp(embeddings_matrix, scale. = TRUE)
    coords <- pca$x[,1:2]
    colnames(coords) <- c("Dim1", "Dim2")
      coords
  } else if (method == "umap") {
      log_info("Computing UMAP")
      umap_result <- umap(embeddings_matrix, n_threads = future::nbrOfWorkers())
      coords <- umap_result$layout
    colnames(coords) <- c("Dim1", "Dim2")
      coords
  } else if (method == "tsne") {
      log_info("Computing t-SNE")
      tsne_result <- Rtsne(embeddings_matrix, perplexity = 30, num.threads = future::nbrOfWorkers())
      coords <- tsne_result$Y
    colnames(coords) <- c("Dim1", "Dim2")
      coords
    }
  }, error = function(e) {
    log_error(sprintf("Error in %s reduction: %s", method, e$message))
    stop(sprintf("Failed to compute %s reduction", method))
  })
  
  log_info(sprintf("Completed %s reduction", method))
  return(as.data.frame(result))
}

# Function to calculate centroid of a group of points
calculate_centroid <- function(points) {
  centroid <- colMeans(points[, c("Dim1", "Dim2")])
  return(centroid)
}

# Function to calculate distances to centroid
calculate_distances_to_centroid <- function(points, centroid) {
  sqrt(rowSums(sweep(points, 2, centroid)^2))
}

# Function to create density plot for centroid distances
create_centroid_density_plot <- function(distances, rounds, reference_type) {
  plot_data <- data.frame(
    distance = distances,
    round = as.factor(rounds)
  )
  
  # Check if we have multiple rounds for comparison
  unique_rounds <- unique(rounds)
  if (length(unique_rounds) > 1) {
    # Perform KS test only if we have multiple rounds
    ks_test <- ks.test(
      distances[rounds == min(unique_rounds)],
      distances[rounds == max(unique_rounds)]
    )
    subtitle <- sprintf("KS test p-value: %.3g", ks_test$p.value)
  } else {
    subtitle <- sprintf("Only Round %s data available", unique_rounds)
  }
  
  ggplot(plot_data, aes(x = distance, fill = round)) +
    geom_density(alpha = 0.4) +
    geom_density(aes(color = round), fill = NA, linewidth = 0.8) +
    scale_fill_manual(values = round_colors) +
    scale_color_manual(values = round_colors) +
    labs(
      title = sprintf("Distance to %s centroid", reference_type),
      subtitle = subtitle,
      x = "Distance to centroid",
      y = "Density",
      fill = "Round",
      color = "Round"
    ) +
    adaptyv_theme() +
    guides(color = "none")
}

# Function to create embedding plot
create_embedding_plot <- function(plot_data, method, ref_sequences = NULL) {
  # Get available rounds
  available_rounds <- unique(plot_data$round)
  round_subtitle <- if (length(available_rounds) > 1) {
    "Colored by round"
  } else {
    sprintf("Round %s data", available_rounds)
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2)) +
    # Add points
    geom_point(aes(color = as.factor(round)), alpha = 0.6, size = 2)
  
  # Add convex hulls for reference-like sequences if provided
  if (!is.null(ref_sequences)) {
    # Calculate similarities for each reference type
    for (ref_type in names(ref_sequences)) {
      similarities <- sapply(plot_data$sequence, function(seq) {
        calculate_sequence_similarity(seq, reference_sequences[[ref_type]])
      })
      
      # Use consistent threshold of 0.5 for selecting similar sequences
      threshold <- 0.8
      similar_sequences <- plot_data$sequence[similarities >= threshold]
      
      # Only create hull if we have enough similar sequences
      if (length(similar_sequences) > 2) {
        hull_data <- plot_data[plot_data$sequence %in% similar_sequences, ]
        hull_color <- switch(ref_type,
                           "egf" = "#FFB6C1",      # Light pink
                           "tgf" = "#B0E0E6",      # Light blue
                           "cetuximab" = "#FFE4B5"  # Light orange/peach
        )
        # Calculate convex hull
        hull_points <- hull_data[chull(hull_data$Dim1, hull_data$Dim2), ]
        # Add hull as polygon
        p <- p + geom_polygon(data = hull_points,
                             aes(x = Dim1, y = Dim2),
                             fill = hull_color,
                             alpha = 0.2,
                             color = hull_color,
                             linewidth = 0.5)
        # Add label at centroid
        centroid <- data.frame(
          x = mean(hull_points$Dim1),
          y = mean(hull_points$Dim2),
          label = toupper(ref_type)
        )
        p <- p + geom_text(data = centroid,
                          aes(x = x, y = y, label = label),
                          color = hull_color,
                          size = 3,
                          fontface = "bold")
      }
    }
  }
  
  # Add styling
  p <- p +
    scale_color_manual(values = round_colors, name = "Round") +
    labs(
      title = sprintf("%s projection of protein sequences", toupper(method)),
      subtitle = round_subtitle,
      x = sprintf("%s dimension 1", toupper(method)),
      y = sprintf("%s dimension 2", toupper(method))
    ) +
    adaptyv_theme()
  
  return(p)
}

# Function to calculate distances to reference centroid
calculate_reference_distances <- function(embeddings, reference_sequences, metadata) {
  # Function to get centroid for a reference type using sequences from both rounds
  get_reference_centroid <- function(ref_type) {
    # Get sequences similar to reference
    similarities <- sapply(metadata$sequence, function(seq) {
      calculate_sequence_similarity(seq, reference_sequences[[ref_type]])
    })
    
    # Use consistent threshold of 0.8
    threshold <- 0.8
    
    # Get sequences from each round
    round1_sequences <- metadata$sequence[similarities >= threshold & metadata$round == 1]
    round2_sequences <- metadata$sequence[similarities >= threshold & metadata$round == 2]
    
    # If not enough sequences in either round, get top sequences from each round
    if (length(round1_sequences) < 3 || length(round2_sequences) < 3) {
      log_warn(sprintf("Not enough sequences with similarity >= %.1f for %s-like centroid in both rounds. Using top sequences from each round.", 
                      threshold, ref_type))
      
      # Get top sequences from each round
      round1_similarities <- similarities[metadata$round == 1]
      round2_similarities <- similarities[metadata$round == 2]
      
      round1_sequences <- metadata$sequence[metadata$round == 1][order(round1_similarities, decreasing = TRUE)][1:3]
      round2_sequences <- metadata$sequence[metadata$round == 2][order(round2_similarities, decreasing = TRUE)][1:3]
    }
    
    ref_like_sequences <- c(round1_sequences, round2_sequences)
    
    # Get embeddings for reference-like sequences
    ref_embeddings <- embeddings[embeddings$sequence %in% ref_like_sequences, -1]
    
    # Log number of sequences used for centroid calculation from each round
    log_info(sprintf("Using %d sequences from Round 1 and %d sequences from Round 2 for %s-like centroid", 
                    sum(metadata$sequence[metadata$round == 1] %in% ref_like_sequences),
                    sum(metadata$sequence[metadata$round == 2] %in% ref_like_sequences),
                    ref_type))
    
    # Compute centroid
    colMeans(as.matrix(ref_embeddings))
  }
  
  # Calculate centroids for each reference type
  centroids <- list(
    egf = get_reference_centroid("egf"),
    tgf = get_reference_centroid("tgf"),
    cetuximab = get_reference_centroid("cetuximab")
  )
  
  # Calculate distances to each centroid
  embedding_matrix <- as.matrix(embeddings[, -1])  # Exclude sequence column
  distances <- list()
  
  for (ref_type in names(centroids)) {
    distances[[ref_type]] <- sqrt(rowSums(sweep(embedding_matrix, 2, centroids[[ref_type]])^2))
  }
  
  # Combine distances with metadata
  result <- data.frame(
    sequence = embeddings$sequence,
    egf_distance = distances$egf,
    tgf_distance = distances$tgf,
    antibody_distance = distances$cetuximab
  )
  
  merge(result, metadata[, c("sequence", "round")], by = "sequence")
}

# Function to create reference distance distribution plots
create_reference_distance_plots <- function(distances_df, output_dir) {
  ref_types <- c("egf", "tgf", "antibody")
  
  for (ref_type in ref_types) {
    dist_col <- paste0(ref_type, "_distance")
    
    # Create density plot
    plot_data <- data.frame(
      distance = distances_df[[dist_col]],
      round = as.factor(distances_df$round)
    )
    
    # Check if we have data for both rounds
    round1_data <- distances_df[[dist_col]][distances_df$round == 1]
    round2_data <- distances_df[[dist_col]][distances_df$round == 2]
    
    # Only perform KS test if we have data for both rounds
    subtitle <- if (length(round1_data) > 0 && length(round2_data) > 0) {
      ks_test <- ks.test(round1_data, round2_data)
      sprintf("KS test p-value: %.3g", ks_test$p.value)
    } else {
      "Insufficient data for statistical comparison"
    }
    
    p <- ggplot(plot_data, aes(x = distance, fill = round)) +
      geom_density(alpha = 0.4) +
      geom_density(aes(color = round), fill = NA, linewidth = 0.8) +
      scale_fill_manual(values = round_colors) +
      scale_color_manual(values = round_colors) +
      labs(
        title = sprintf("Distance to %s-like centroid", toupper(ref_type)),
        subtitle = subtitle,
        x = "Distance in embedding space",
        y = "Density",
        fill = "Round",
        color = "Round"
      ) +
      adaptyv_theme() +
      guides(color = "none")
    
    ggsave(
      file.path(output_dir, sprintf("reference_distance_%s.svg", ref_type)),
      p, width = 10, height = 8
    )
  }
}

main <- function() {
  tic("Total execution time")
  
  # Read metadata first to get filtered sequences
  tic("Reading metadata")
  tryCatch({
    metadata <- fread(opt$metadata)
    
    # Ensure column names are consistent
    setnames(metadata, tolower(names(metadata)))
    
    # Convert round to numeric if it isn't already
    metadata[, round := as.numeric(as.character(round))]
    
    # Log total sequences per round in metadata
    log_info(sprintf("Total sequences in metadata - Round 1: %d, Round 2: %d", 
                    sum(metadata$round == 1),
                    sum(metadata$round == 2)))
    
    # Filter by round if specified
    if (opt$round != "both") {
      round_num <- as.numeric(opt$round)
      metadata <- metadata[round == round_num]
      log_info(sprintf('Filtered data for round %d', round_num))
    }

    # Filter for tested designs if specified
    if (opt$tested_only) {
      log_info('Filtering for tested designs only...')
      metadata <- metadata[!is.na(selected) & selected %in% c("Adaptyv selection", "Top 100")]
      log_info(sprintf("Tested designs - Round 1: %d, Round 2: %d", 
                      sum(metadata$round == 1, na.rm = TRUE),
                      sum(metadata$round == 2, na.rm = TRUE)))
    } else {
      log_info('Including all sequences...')
      log_info(sprintf("Total sequences - Round 1: %d, Round 2: %d", 
                      sum(metadata$round == 1, na.rm = TRUE),
                      sum(metadata$round == 2, na.rm = TRUE)))
    }
  }, error = function(e) {
    log_error(sprintf("Error reading metadata: %s", e$message))
    stop("Failed to process metadata")
  })
  toc()
  
  # Read embeddings and compute reference distances
  tic("Computing reference distances")
  tryCatch({
    # Read embeddings
    embeddings <- read_embeddings(opt$embeddings, metadata)
    
    # Calculate distances to reference centroids and get reference-like sequences
    log_info("Calculating distances to reference centroids")
    ref_sequences <- list()
    for (ref_type in c("egf", "tgf", "cetuximab")) {
      similarities <- sapply(metadata$sequence, function(seq) {
        calculate_sequence_similarity(seq, reference_sequences[[ref_type]])
      })
      threshold <- 0.8  # Use 0.5 for all reference types
      ref_sequences[[ref_type]] <- metadata$sequence[similarities >= threshold]
      log_info(sprintf("Found %d sequences with similarity >= %.1f for %s-like centroid",
                      length(ref_sequences[[ref_type]]), threshold, ref_type))
    }
    
    # Perform dimensionality reduction
    coords <- reduce_dimensions(embeddings[, -1], opt$method)
    
    # Ensure all sequences have round information
    sequence_rounds <- metadata[, .(sequence, round)]
    plot_data <- data.frame(
      sequence = embeddings$sequence,
      coords
    )
    plot_data <- merge(plot_data, sequence_rounds, by = "sequence", all.x = TRUE)
    
    # Log any sequences without round information
    missing_rounds <- sum(is.na(plot_data$round))
    if (missing_rounds > 0) {
      log_warn(sprintf("%d sequences are missing round information", missing_rounds))
    }
    
    # Create embedding plot with convex hulls
    p <- create_embedding_plot(plot_data, opt$method, ref_sequences)
    
    # Save embedding plot
    dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
    ggsave(
      file.path(opt$output, sprintf("embedding_%s.svg", opt$method)),
      p, width = 10, height = 8
    )
    
    # Calculate and plot reference distances
    distances_df <- calculate_reference_distances(embeddings, reference_sequences, metadata)
    create_reference_distance_plots(distances_df, opt$output)
    
  }, error = function(e) {
    log_error(sprintf("Error in reference distance calculation: %s", e$message))
    stop("Failed to compute reference distances")
  })
  toc()
  
  toc()  # Total execution time
  log_info("All plots have been generated successfully")
}

# Helper function to read embeddings
read_embeddings <- function(embeddings_path, metadata) {
  if (grepl("\\.h5$", embeddings_path)) {
    log_info("Reading H5 file")
    h5_datasets <- h5ls(embeddings_path)
    dataset_names <- h5_datasets$name
    
    # First, create a mapping between dataset names and sequences
    sequence_mapping <- data.frame(
      dataset = dataset_names,
      sequence = metadata$sequence[1:length(dataset_names)]
    )
    
    # Log sequence mapping stats
    log_info(sprintf("Found %d sequences in h5 file", length(dataset_names)))
    log_info(sprintf("Matched with %d sequences from metadata", sum(sequence_mapping$sequence %in% metadata$sequence)))
    
    # Get embedding dimension from first sequence
    first_embedding <- h5read(embeddings_path, dataset_names[1])
    embedding_dim <- nrow(first_embedding)  # Use nrow instead of ncol since it's a 2560x1 matrix
    log_info(sprintf("Embedding dimension: %d", embedding_dim))
    
    embeddings_matrix <- matrix(0, nrow = length(dataset_names), ncol = embedding_dim)
    chunk_size <- ceiling(length(dataset_names) / future::nbrOfWorkers())
    chunks <- split(seq_along(dataset_names), 
                   ceiling(seq_along(dataset_names) / chunk_size))
    
    log_info("Reading embeddings in parallel chunks")
    embeddings_list <- future_map(chunks, function(chunk_indices) {
      chunk_matrix <- matrix(0, nrow = length(chunk_indices), ncol = embedding_dim)
      for (i in seq_along(chunk_indices)) {
        idx <- chunk_indices[i]
        embedding <- h5read(embeddings_path, dataset_names[idx])
        chunk_matrix[i,] <- as.numeric(t(embedding))  # Transpose the 2560x1 matrix to 1x2560
      }
      return(chunk_matrix)
    }, .progress = TRUE)
    
    embeddings_matrix <- do.call(rbind, embeddings_list)
    log_info("Combined all embedding chunks")
    
    # Only include sequences that are in metadata
    valid_sequences <- sequence_mapping$sequence %in% metadata$sequence
    result <- data.frame(
      sequence = sequence_mapping$sequence[valid_sequences],
      as.data.frame(embeddings_matrix[valid_sequences,])
    )
    
    log_info(sprintf("Returning embeddings for %d sequences that match metadata", nrow(result)))
    return(result)
  } else if (grepl("\\.parquet$", embeddings_path)) {
    embeddings <- arrow::read_parquet(embeddings_path)
    if (is.list(embeddings[[2]])) {
      embedding_matrix <- do.call(rbind, embeddings[[2]])
      data.frame(
        sequence = embeddings[[1]],
        as.data.frame(embedding_matrix)
      )
    } else {
      embeddings
    }
  } else {
    fread(embeddings_path)
  }
}

if (sys.nframe() == 0) {
  main()
}
