## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE, 
  warning = FALSE, 
  message = FALSE,
  fig.width = 8, 
  fig.height = 6,
  collapse = TRUE,
  comment = "#>",
  cache = FALSE
)
library(topolow)

# Function to check and load packages with informative messages
check_and_load_package <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    message("Package ", package_name, " not found. Please install it with:")
    message("install.packages('", package_name, "')")
    return(FALSE)
  } else {
    library(package_name, character.only = TRUE)
    return(TRUE)
  }
}

# Load required libraries
required_packages <- c("topolow", "ggplot2", "dplyr", "reshape2", 
                      "scales", "MASS")

# Check for optional packages
optional_packages <- c("smacof", "RANN")

missing_packages <- character(0)
for(pkg in required_packages) {
  if(!check_and_load_package(pkg)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

# Load optional packages quietly
for(pkg in optional_packages) {
  check_and_load_package(pkg)
}

if(length(missing_packages) > 0) {
  stop("Please install missing packages: ", paste(missing_packages, collapse = ", "))
}

# Set seed for reproducibility
set.seed(42)

## ----utility_functions, eval=FALSE--------------------------------------------
# # Enhanced statistical testing with effect sizes
# enhanced_statistical_comparison <- function(topolow_errors, other_errors, method_name) {
#   valid_topolow <- topolow_errors[!is.na(topolow_errors)]
#   valid_other <- other_errors[!is.na(other_errors)]
# 
#   if(length(valid_other) < 2 || length(valid_topolow) < 2) {
#     cat("Insufficient data for statistical comparison with", method_name, "\n")
#     return(NULL)
#   }
# 
#   # T-test
#   t_result <- t.test(valid_topolow, valid_other)
# 
#   # Effect size (Cohen's d)
#   pooled_sd <- sqrt(((length(valid_topolow) - 1) * var(valid_topolow) +
#                      (length(valid_other) - 1) * var(valid_other)) /
#                     (length(valid_topolow) + length(valid_other) - 2))
# 
#   cohens_d <- (mean(valid_topolow) - mean(valid_other)) / pooled_sd
# 
#   cat("\n--- Statistical Comparison (Topolow vs", method_name, ") ---\n")
#   cat("- Welch's t-statistic:", round(t_result$statistic, 3), "\n")
#   cat("- p-value:", format(t_result$p.value, scientific = TRUE, digits = 3), "\n")
#   cat("- Cohen's d (effect size):", round(cohens_d, 3), "\n")
#   cat("- Effect size interpretation:",
#       if(abs(cohens_d) < 0.2) "Negligible"
#       else if(abs(cohens_d) < 0.5) "Small"
#       else if(abs(cohens_d) < 0.8) "Medium"
#       else "Large", "\n")
# 
#   return(list(t_test = t_result, cohens_d = cohens_d))
# }
# 
# # Calculate data quality metrics
# calculate_data_quality_metrics <- function(distance_matrix) {
#   total_possible <- nrow(distance_matrix) * (nrow(distance_matrix) - 1) / 2
#   total_available <- sum(!is.na(distance_matrix[upper.tri(distance_matrix)]))
#   completeness <- total_available / total_possible
# 
#   available_distances <- distance_matrix[!is.na(distance_matrix)]
# 
#   metrics <- list(
#     completeness = completeness,
#     n_objects = nrow(distance_matrix),
#     n_available_distances = total_available,
#     distance_range = range(available_distances),
#     distance_mean = mean(available_distances),
#     distance_sd = sd(available_distances),
#     distance_cv = sd(available_distances) / mean(available_distances)
#   )
#   return(metrics)
# }
# 
# # Coordinate validation
# validate_coordinates <- function(coords, method_name, n_objects, target_dims) {
#   if(is.null(coords)) {
#     cat("WARNING:", method_name, "returned NULL coordinates\n")
#     return(FALSE)
#   }
#   if(any(!is.finite(coords))) {
#     cat("WARNING:", method_name, "has non-finite coordinates\n")
#     return(FALSE)
#   }
#   if(nrow(coords) != n_objects || ncol(coords) != target_dims) {
#     cat("WARNING:", method_name, "has incorrect dimensions\n")
#     return(FALSE)
#   }
#   return(TRUE)
# }
# 
# # Missing data imputation
# improved_missing_data_imputation <- function(dist_matrix) {
#   available_distances <- dist_matrix[!is.na(dist_matrix)]
#   if(length(available_distances) == 0) {
#     median_distance <- 1.0
#   } else {
#     median_distance <- median(available_distances, na.rm = TRUE)
#   }
#   dist_matrix[is.na(dist_matrix)] <- median_distance
#   return(list(matrix = dist_matrix, imputation_value = median_distance))
# }

## ----distorted_euclidean_generator, eval=FALSE--------------------------------
# #' Generate Non-Euclidean Data by Distorting Clustered Points
# #'
# #' Creates non-Euclidean data through systematic distortion of clustered high-dimensional points.
# #' This method is particularly effective for testing robustness to violations of the triangle inequality.
# #'
# #' @param n_objects Number of objects to generate
# #' @param initial_dims Dimensionality of the initial latent space
# #' @param n_clusters Number of clusters to form
# #' @param noise_factor Magnitude of asymmetric noise to add
# #' @param missing_fraction Proportion of distances to set to NA
# #' @return List containing complete and incomplete distance matrices
# generate_distorted_euclidean_data <- function(n_objects = 50,
#                                               initial_dims = 10,
#                                               n_clusters = 8,
#                                               noise_factor = 0.3,
#                                               missing_fraction = 0.30) {
# 
#   object_names <- paste0("Object_", sprintf("%02d", 1:n_objects))
# 
#   # Generate structured high-dimensional coordinates
#   cluster_size <- n_objects %/% n_clusters
#   cluster_centers <- matrix(rnorm(n_clusters * initial_dims, mean = 0, sd = 4),
#                            nrow = n_clusters, ncol = initial_dims)
#   initial_coords <- matrix(0, nrow = n_objects, ncol = initial_dims)
#   rownames(initial_coords) <- object_names
# 
#   for(i in 1:n_objects) {
#     cluster_id <- ((i - 1) %/% cluster_size) + 1
#     if(cluster_id > n_clusters) cluster_id <- n_clusters
#     initial_coords[i, ] <- cluster_centers[cluster_id, ] +
#                           rnorm(initial_dims, mean = 0, sd = 1.5)
#   }
# 
#   # Calculate foundational Euclidean distances
#   euclidean_distances <- as.matrix(dist(initial_coords, method = "euclidean"))
# 
#   # Apply non-linear transformations
#   dist_quantiles <- quantile(euclidean_distances[upper.tri(euclidean_distances)],
#                             c(0.33, 0.66))
#   transform_1 <- euclidean_distances
#   transform_1[euclidean_distances <= dist_quantiles[1]] <-
#     euclidean_distances[euclidean_distances <= dist_quantiles[1]]^1.3
#   transform_1[euclidean_distances > dist_quantiles[1] &
#              euclidean_distances <= dist_quantiles[2]] <-
#     euclidean_distances[euclidean_distances > dist_quantiles[1] &
#                        euclidean_distances <= dist_quantiles[2]]^1.6
#   transform_1[euclidean_distances > dist_quantiles[2]] <-
#     euclidean_distances[euclidean_distances > dist_quantiles[2]]^1.8
# 
#   # Add asymmetric noise
#   asymmetric_noise <- transform_1 * noise_factor *
#                      matrix(runif(n_objects^2, -1, 1), nrow = n_objects)
#   asymmetric_noise <- (asymmetric_noise + t(asymmetric_noise)) / 2
#   transform_2 <- transform_1 + asymmetric_noise
#   transform_2[transform_2 < 0] <- 0.01
# 
#   # Create final non-Euclidean distance matrix
#   complete_non_euclidean_distances <- transform_2
#   diag(complete_non_euclidean_distances) <- 0
# 
#   # Introduce missing values
#   total_unique_pairs <- n_objects * (n_objects - 1) / 2
#   target_missing_pairs <- round(total_unique_pairs * missing_fraction)
# 
#   upper_tri_indices <- which(upper.tri(complete_non_euclidean_distances), arr.ind = TRUE)
#   missing_pair_indices <- sample(nrow(upper_tri_indices), target_missing_pairs)
# 
#   incomplete_distances <- complete_non_euclidean_distances
#   incomplete_distances[upper_tri_indices[missing_pair_indices,]] <- NA
#   incomplete_distances[upper_tri_indices[missing_pair_indices, c(2,1)]] <- NA
# 
#   return(list(
#     complete_distances = complete_non_euclidean_distances,
#     incomplete_distances = incomplete_distances,
#     object_names = object_names,
#     method = "Synthesized non-Euclidean"
#   ))
# }

## ----swiss_roll_generator, eval=FALSE-----------------------------------------
# #' Generate Non-Euclidean Data from a Swiss Roll Manifold
# #'
# #' Creates data points on a Swiss Roll manifold where geodesic distances
# #' are inherently non-Euclidean. Includes "tunneling" effects common in real-world manifold data.
# #'
# #' @param n_objects Number of points on the manifold
# #' @param noise Standard deviation of Gaussian noise added to points
# #' @param tunnel_fraction Fraction of distances to replace with Euclidean shortcuts
# #' @param missing_fraction Proportion of final distances to set to NA
# #' @return List containing complete and incomplete distance matrices
# generate_swiss_roll_data <- function(n_objects = 50,
#                                      noise = 0.05,
#                                      tunnel_fraction = 0.05,
#                                      missing_fraction = 0.30) {
#   # Check if required packages are available
#   if(!requireNamespace("RANN", quietly = TRUE)) {
#     stop("RANN package is required for this function. Please install it.")
#   }
#   if(!requireNamespace("igraph", quietly = TRUE)) {
#     stop("igraph package is required for this function. Please install it.")
#   }
# 
#   # Generate points on the Swiss Roll
#   t <- 1.5 * pi * (1 + 2 * runif(n_objects))
#   height <- 21 * runif(n_objects)
# 
#   x <- t * cos(t)
#   y <- height
#   z <- t * sin(t)
# 
#   points_3d <- data.frame(x = x, y = y, z = z) + rnorm(n_objects * 3, sd = noise)
#   object_names <- paste0("Object_", sprintf("%02d", 1:n_objects))
#   rownames(points_3d) <- object_names
# 
#   # Calculate geodesic distances via k-NN graph
#   if(requireNamespace("RANN", quietly = TRUE)) {
#     knn_graph <- RANN::nn2(points_3d, k = min(10, n_objects-1))$nn.idx
#     adj_matrix <- matrix(0, n_objects, n_objects)
#     for (i in 1:n_objects) {
#       for (j_idx in 2:min(10, n_objects)) {
#         if(j_idx <= ncol(knn_graph)) {
#           j <- knn_graph[i, j_idx]
#           dist_ij <- dist(rbind(points_3d[i,], points_3d[j,]))
#           adj_matrix[i, j] <- dist_ij
#           adj_matrix[j, i] <- dist_ij
#         }
#       }
#     }
# 
#     g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
#     geodesic_distances <- igraph::distances(g)
#     geodesic_distances[is.infinite(geodesic_distances)] <-
#       max(geodesic_distances[is.finite(geodesic_distances)]) * 1.5
#   } else {
#     # Fallback to Euclidean if RANN not available
#     geodesic_distances <- as.matrix(dist(points_3d))
#     warning("RANN package not available. Using Euclidean distances as approximation.")
#   }
# 
#   # Introduce "tunnels" or "short-circuits"
#   euclidean_distances <- as.matrix(dist(points_3d))
#   n_tunnels <- round(tunnel_fraction * n_objects * (n_objects - 1) / 2)
# 
#   upper_tri_indices <- which(upper.tri(geodesic_distances), arr.ind = TRUE)
#   tunnel_indices <- sample(nrow(upper_tri_indices), min(n_tunnels, nrow(upper_tri_indices)))
# 
#   complete_distances <- geodesic_distances
#   if(length(tunnel_indices) > 0) {
#     for (k in tunnel_indices) {
#       i <- upper_tri_indices[k, 1]
#       j <- upper_tri_indices[k, 2]
#       complete_distances[i, j] <- complete_distances[j, i] <- euclidean_distances[i, j]
#     }
#   }
#   diag(complete_distances) <- 0
# 
#   # Introduce missing values
#   total_unique_pairs <- n_objects * (n_objects - 1) / 2
#   target_missing_pairs <- round(total_unique_pairs * missing_fraction)
# 
#   missing_pair_indices <- sample(nrow(upper_tri_indices),
#                                 min(target_missing_pairs, nrow(upper_tri_indices)))
# 
#   incomplete_distances <- complete_distances
#   if(length(missing_pair_indices) > 0) {
#     incomplete_distances[upper_tri_indices[missing_pair_indices,]] <- NA
#     incomplete_distances[upper_tri_indices[missing_pair_indices, c(2,1)]] <- NA
#   }
# 
#   rownames(complete_distances) <- object_names
#   colnames(complete_distances) <- object_names
#   rownames(incomplete_distances) <- object_names
#   colnames(incomplete_distances) <- object_names
# 
#   return(list(
#     complete_distances = complete_distances,
#     incomplete_distances = incomplete_distances,
#     object_names = object_names,
#     method = "Swiss Roll Manifold"
#   ))
# }

## ----analysis_pipeline, eval=FALSE--------------------------------------------
# #' Comprehensive analysis pipeline for embedding method comparison
# #'
# #' @param data_gen_func Function to generate the dataset
# #' @param dataset_name Name for the dataset (for labeling)
# #' @param n_objects Number of objects in the dataset
# #' @param n_runs Number of runs for stochastic methods
# run_comprehensive_analysis <- function(data_gen_func, dataset_name, n_objects = 50, n_runs = 50) {
# 
#   cat("=== ANALYSIS FOR", toupper(dataset_name), "===\n")
# 
#   # Step 1: Data Generation
#   cat("\n1. Generating", dataset_name, "dataset...\n")
#   data_gen_output <- data_gen_func(n_objects = n_objects, missing_fraction = 0.3)
# 
#   complete_distances_for_evaluation <- data_gen_output$complete_distances
#   incomplete_distances_for_embedding <- data_gen_output$incomplete_distances
#   object_names <- data_gen_output$object_names
# 
#   actual_missing_percentage <- sum(is.na(incomplete_distances_for_embedding)) /
#                               (n_objects * (n_objects-1)) * 100
# 
#   data_quality <- calculate_data_quality_metrics(incomplete_distances_for_embedding)
# 
#   cat("Generated dataset with", n_objects, "objects and",
#       round(actual_missing_percentage, 1), "% missing values\n")
# 
#   # Step 2: Assess Non-Euclidean Character
#   cat("\n2. Assessing non-Euclidean character...\n")
# 
#   D_squared <- complete_distances_for_evaluation^2
#   n <- nrow(D_squared)
#   J <- diag(n) - (1/n) * matrix(1, n, n)
#   B <- -0.5 * J %*% D_squared %*% J
#   eigenvals <- eigen(B, symmetric = TRUE)$values
# 
#   numerical_tolerance <- 1e-12
#   positive_eigenvals <- eigenvals[eigenvals > numerical_tolerance]
#   negative_eigenvals <- eigenvals[eigenvals < -numerical_tolerance]
# 
#   total_variance <- sum(abs(eigenvals))
#   negative_variance <- sum(abs(negative_eigenvals))
#   positive_variance <- sum(positive_eigenvals)
# 
#   if(positive_variance > 0) {
#     deviation_score <- negative_variance / positive_variance
#   } else {
#     deviation_score <- 1.0
#   }
# 
#   negative_variance_fraction <- negative_variance / total_variance
# 
#   cumulative_positive_variance <- cumsum(positive_eigenvals) / positive_variance
#   dims_for_90_percent <- which(cumulative_positive_variance >= 0.90)[1]
#   if(is.na(dims_for_90_percent)) dims_for_90_percent <- length(positive_eigenvals)
# 
#   cat("Non-Euclidean deviation score:", round(deviation_score, 4), "\n")
#   cat("Dimensions for 90% variance:", dims_for_90_percent, "\n")
# 
#   # Step 3: Parameter Optimization using Euclidify
#   cat("\n3. Running automated parameter optimization...\n")
# 
#   # Create temporary directory for optimization
#   temp_dir <- tempfile()
#   dir.create(temp_dir, recursive = TRUE)
# 
#   euclidify_results <- tryCatch({
#     Euclidify(
#       dissimilarity_matrix = incomplete_distances_for_embedding,
#       output_dir = temp_dir,
#       ndim_range = c(2, min(10, dims_for_90_percent + 3)),
#       n_initial_samples = 20,  # Reduced for vignette speed
#       n_adaptive_samples = 40,  # Reduced for vignette speed
#       folds = 20,
#       mapping_max_iter = 300,  # Reduced for speed
#       clean_intermediate = TRUE,
#       verbose = "off",
#       fallback_to_defaults = TRUE,
#       save_results = FALSE
#     )
#   }, error = function(e) {
#     cat("Euclidify failed, using default parameters\n")
#     NULL
#   })
# 
#   # Clean up temp directory
#   unlink(temp_dir, recursive = TRUE)
# 
#   if(!is.null(euclidify_results)) {
#     optimal_params <- euclidify_results$optimal_params
#     euclidify_positions <- euclidify_results$positions
#     target_dims <- optimal_params$ndim
#     cat("Optimal parameters found - dimensions:", target_dims, "\n")
#   } else {
#     # Fallback parameters
#     target_dims <- max(2, min(5, dims_for_90_percent))
#     optimal_params <- list(
#       ndim = target_dims,
#       k0 = 5.0,
#       cooling_rate = 0.01,
#       c_repulsion = 0.02
#     )
#     euclidify_positions <- NULL
#     cat("Using fallback parameters - dimensions:", target_dims, "\n")
#   }
# 
#   # Step 4: Three-Method Comparison
#   cat("\n4. Running three-method comparison...\n")
# 
#   # Initialize storage
#   topolow_results <- vector("list", n_runs)
#   iterative_mds_results <- vector("list", n_runs)
#   classical_mds_result <- NULL
# 
#   topolow_errors <- numeric(n_runs)
#   iterative_mds_errors <- numeric(n_runs)
#   classical_mds_error <- NA
# 
#   # Topolow runs
#   cat("Running Topolow embeddings...\n")
#   for(i in 1:n_runs) {
#     if(i == 1 && !is.null(euclidify_positions)) {
#       # Use Euclidify result for first run
#       topolow_coords <- euclidify_positions
#     } else {
#       # Run new embedding
#       topolow_result <- tryCatch({
#         euclidean_embedding(
#           dissimilarity_matrix = incomplete_distances_for_embedding,
#           ndim = optimal_params$ndim,
#           mapping_max_iter = 300,
#           k0 = optimal_params$k0,
#           cooling_rate = optimal_params$cooling_rate,
#           c_repulsion = optimal_params$c_repulsion,
#           relative_epsilon = 1e-6,
#           convergence_counter = 3,
#           verbose = FALSE
#         )
#       }, error = function(e) NULL)
# 
#       if(!is.null(topolow_result)) {
#         topolow_coords <- topolow_result$positions
#       } else {
#         topolow_coords <- NULL
#       }
#     }
# 
#     if(!is.null(topolow_coords)) {
#       topolow_coords <- topolow_coords[order(row.names(topolow_coords)), ]
# 
#       if(validate_coordinates(topolow_coords, "Topolow", n_objects, target_dims)) {
#         topolow_coords <- scale(topolow_coords, center = TRUE, scale = FALSE)
#         topolow_results[[i]] <- topolow_coords
# 
#         embedded_distances <- as.matrix(dist(topolow_coords))
#         rownames(embedded_distances) <- rownames(topolow_coords)
#         colnames(embedded_distances) <- rownames(topolow_coords)
# 
#         valid_mask <- !is.na(complete_distances_for_evaluation)
#         distance_errors <- abs(complete_distances_for_evaluation[valid_mask] -
#                              embedded_distances[valid_mask])
#         topolow_errors[i] <- sum(distance_errors)
#       } else {
#         topolow_results[[i]] <- NULL
#         topolow_errors[i] <- NA
#       }
#     } else {
#       topolow_results[[i]] <- NULL
#       topolow_errors[i] <- NA
#     }
#   }
# 
#   # Classical MDS
#   cat("Running Classical MDS...\n")
#   imputation_result <- improved_missing_data_imputation(incomplete_distances_for_embedding)
#   classical_mds_matrix <- imputation_result$matrix
# 
#   tryCatch({
#     classical_mds_result_raw <- cmdscale(classical_mds_matrix, k = target_dims, eig = TRUE)
#     classical_mds_coords <- classical_mds_result_raw$points
#     rownames(classical_mds_coords) <- object_names
#     classical_mds_coords <- classical_mds_coords[order(row.names(classical_mds_coords)), ]
# 
#     if(validate_coordinates(classical_mds_coords, "Classical MDS", n_objects, target_dims)) {
#       classical_mds_coords <- scale(classical_mds_coords, center = TRUE, scale = FALSE)
#       classical_mds_result <- classical_mds_coords
# 
#       embedded_distances <- as.matrix(dist(classical_mds_coords))
#       rownames(embedded_distances) <- rownames(classical_mds_coords)
#       colnames(embedded_distances) <- rownames(classical_mds_coords)
# 
#       valid_mask <- !is.na(complete_distances_for_evaluation)
#       distance_errors <- abs(complete_distances_for_evaluation[valid_mask] -
#                            embedded_distances[valid_mask])
#       classical_mds_error <- sum(distance_errors)
#     }
#   }, error = function(e) {
#     cat("Classical MDS failed:", e$message, "\n")
#   })
# 
#   # Iterative MDS
#   cat("Running Iterative MDS...\n")
#   for(i in 1:n_runs) {
#     tryCatch({
#       if(requireNamespace("smacof", quietly = TRUE)) {
#         max_dist <- max(classical_mds_matrix, na.rm = TRUE)
#         scaled_matrix <- classical_mds_matrix / max_dist
# 
#         iterative_mds_result_raw <- smacof::smacofSym(
#           delta = scaled_matrix,
#           ndim = target_dims,
#           type = "interval",
#           init = "torgerson",
#           verbose = FALSE,
#           itmax = 1000,
#           eps = 1e-5
#         )
# 
#         iterative_mds_coords <- iterative_mds_result_raw$conf
#         current_max_dist <- max(dist(iterative_mds_coords))
#         if(current_max_dist > 0) {
#           scale_factor <- max_dist / current_max_dist
#           iterative_mds_coords <- iterative_mds_coords * scale_factor
#         }
#       } else {
#         # Fallback to isoMDS
#         iterative_mds_result_raw <- MASS::isoMDS(classical_mds_matrix, k = target_dims, trace = FALSE)
#         iterative_mds_coords <- iterative_mds_result_raw$points
#       }
# 
#       rownames(iterative_mds_coords) <- object_names
#       iterative_mds_coords <- iterative_mds_coords[order(row.names(iterative_mds_coords)), ]
# 
#       if(validate_coordinates(iterative_mds_coords, "Iterative MDS", n_objects, target_dims)) {
#         iterative_mds_coords <- scale(iterative_mds_coords, center = TRUE, scale = FALSE)
#         iterative_mds_results[[i]] <- iterative_mds_coords
# 
#         embedded_distances <- as.matrix(dist(iterative_mds_coords))
#         rownames(embedded_distances) <- rownames(iterative_mds_coords)
#         colnames(embedded_distances) <- rownames(iterative_mds_coords)
# 
#         valid_mask <- !is.na(complete_distances_for_evaluation)
#         distance_errors <- abs(complete_distances_for_evaluation[valid_mask] -
#                              embedded_distances[valid_mask])
#         iterative_mds_errors[i] <- sum(distance_errors)
#       } else {
#         iterative_mds_results[[i]] <- NULL
#         iterative_mds_errors[i] <- NA
#       }
# 
#     }, error = function(e) {
#       iterative_mds_results[[i]] <- NULL
#       iterative_mds_errors[i] <- NA
#     })
#   }
# 
#   # Step 5: Compile Results
#   valid_topolow_results <- !is.na(topolow_errors)
#   valid_iterative_results <- !is.na(iterative_mds_errors)
# 
#   results <- list(
#     dataset_name = dataset_name,
#     data_characteristics = list(
#       n_objects = n_objects,
#       missing_percentage = actual_missing_percentage,
#       deviation_score = deviation_score,
#       dims_90_percent = dims_for_90_percent,
#       negative_variance_fraction = negative_variance_fraction
#     ),
#     optimal_params = optimal_params,
#     topolow_errors = topolow_errors,
#     topolow_results = topolow_results,
#     valid_topolow_results = valid_topolow_results,
#     classical_mds_error = classical_mds_error,
#     classical_mds_result = classical_mds_result,
#     iterative_mds_errors = iterative_mds_errors,
#     iterative_mds_results = iterative_mds_results,
#     valid_iterative_results = valid_iterative_results,
#     complete_distances = complete_distances_for_evaluation,
#     incomplete_distances = incomplete_distances_for_embedding,
#     target_dims = target_dims
#   )
# 
#   return(results)
# }

## ----distorted_analysis, eval=FALSE, results='hide', fig.show='hold'----------
# # Run analysis for Synthesized non-Euclidean data
# distorted_results <- run_comprehensive_analysis(
#   generate_distorted_euclidean_data,
#   "Synthesized non-Euclidean",
#   n_objects = 50,
#   n_runs = 50
# )

## ----distorted_summary, eval=FALSE, echo=FALSE--------------------------------
# cat("=== Synthesized non-Euclidean DATA ANALYSIS SUMMARY ===\n")
# cat("Dataset characteristics:\n")
# cat("- Objects:", distorted_results$data_characteristics$n_objects, "\n")
# cat("- Missing data:", round(distorted_results$data_characteristics$missing_percentage, 1), "%\n")
# cat("- Non-Euclidean deviation score:", round(distorted_results$data_characteristics$deviation_score, 4), "\n")
# cat("- Dimensions for 90% variance:", distorted_results$data_characteristics$dims_90_percent, "\n")
# cat("- Target embedding dimensions:", distorted_results$target_dims, "\n\n")
# 
# cat("Performance Summary:\n")
# if(sum(distorted_results$valid_topolow_results) > 0) {
#   cat("- Topolow mean error:",
#       round(mean(distorted_results$topolow_errors[distorted_results$valid_topolow_results]), 2),
#       "±", round(sd(distorted_results$topolow_errors[distorted_results$valid_topolow_results]), 2), "\n")
# }
# if(!is.na(distorted_results$classical_mds_error)) {
#   cat("- Classical MDS error:", round(distorted_results$classical_mds_error, 2), "\n")
# }
# if(sum(distorted_results$valid_iterative_results) > 0) {
#   cat("- Iterative MDS mean error:", round(mean(distorted_results$iterative_mds_errors[distorted_results$valid_iterative_results]), 2),
#       "±", round(sd(distorted_results$iterative_mds_errors[distorted_results$valid_iterative_results]), 2), "\n")
# }

## ----distorted_eigenvalues, eval=FALSE, echo=FALSE----------------------------
# # Create eigenvalue visualization for distorted data
# D_squared <- distorted_results$complete_distances^2
# n <- nrow(D_squared)
# J <- diag(n) - (1/n) * matrix(1, n, n)
# B <- -0.5 * J %*% D_squared %*% J
# eigenvals <- eigen(B, symmetric = TRUE)$values
# 
# eigenval_df <- data.frame(
#   index = 1:length(eigenvals),
#   eigenvalue = eigenvals,
#   type = ifelse(eigenvals > 1e-12, "Positive",
#                ifelse(eigenvals < -1e-12, "Negative", "Zero"))
# )
# 
# p_eigen_distorted <- ggplot(eigenval_df, aes(x = index, y = eigenvalue, fill = type)) +
#   geom_bar(stat = "identity", alpha = 0.8) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "black") +
#   scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "red", "Zero" = "gray")) +
#   labs(
#     title = "Eigenvalue Spectrum: Synthesized non-Euclidean Data",
#     subtitle = paste("Deviation Score:", round(distorted_results$data_characteristics$deviation_score, 4),
#                      "| Non-Euclidean Variance:", round(distorted_results$data_characteristics$negative_variance_fraction * 100, 1), "%"),
#     x = "Eigenvalue Index (sorted descending)",
#     y = "Eigenvalue",
#     fill = "Eigenvalue Type"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(size = 12, face = "bold"),
#     panel.grid.minor = element_blank()
#   )
# 
# print(p_eigen_distorted)

## ----distorted_eigenvalues-1, echo=FALSE, out.width="100%"--------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/distorted_eigenvalues-1.png")

## ----distorted_performance, eval=FALSE, echo=FALSE----------------------------
# # Create performance comparison for distorted data
# results_df_distorted <- data.frame()
# 
# if(sum(distorted_results$valid_topolow_results) > 0) {
#   topolow_df <- data.frame(
#     Run = which(distorted_results$valid_topolow_results),
#     Method = "Topolow",
#     Total_Error = distorted_results$topolow_errors[distorted_results$valid_topolow_results],
#     Method_Type = "Stochastic",
#     Dataset = "Synthesized non-Euclidean"
#   )
#   results_df_distorted <- rbind(results_df_distorted, topolow_df)
# }
# 
# if(!is.na(distorted_results$classical_mds_error)) {
#   classical_df <- data.frame(
#     Run = 1,
#     Method = "Classical MDS",
#     Total_Error = distorted_results$classical_mds_error,
#     Method_Type = "Deterministic",
#     Dataset = "Synthesized non-Euclidean"
#   )
#   results_df_distorted <- rbind(results_df_distorted, classical_df)
# }
# 
# if(sum(distorted_results$valid_iterative_results) > 0) {
#   iterative_df <- data.frame(
#     Run = which(distorted_results$valid_iterative_results),
#     Method = "Iterative MDS",
#     Total_Error = distorted_results$iterative_mds_errors[distorted_results$valid_iterative_results],
#     Method_Type = "Stochastic",
#     Dataset = "Synthesized non-Euclidean"
#   )
#   results_df_distorted <- rbind(results_df_distorted, iterative_df)
# }
# 
# if(nrow(results_df_distorted) > 0) {
#   p_performance_distorted <- ggplot(results_df_distorted, aes(x = Method, y = Total_Error, fill = Method_Type)) +
#     geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.7, size = 2.5) +
#     scale_fill_manual(values = c("Stochastic" = "steelblue", "Deterministic" = "darkorange")) +
#     labs(
#       title = "Embedding Performance: Synthesized non-Euclidean Data",
#       subtitle = paste("Non-Euclidean Score:", round(distorted_results$data_characteristics$deviation_score, 4),
#                        "| Missing Data:", round(distorted_results$data_characteristics$missing_percentage, 1), "%"),
#       y = "Sum of Absolute Distance Errors (L1 Norm)",
#       x = "Method",
#       fill = "Method Type"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 12, face = "bold"),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
# 
#   print(p_performance_distorted)
# }

## ----distorted_performance-1, echo=FALSE, out.width="100%"--------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/distorted_performance-1.png")

## ----swiss_analysis, eval=FALSE, results='hide', fig.show='hold'--------------
# # Run analysis for Swiss Roll data
# swiss_results <- run_comprehensive_analysis(
#   generate_swiss_roll_data,
#   "Swiss Roll Manifold",
#   n_objects = 50,
#   n_runs = 50
# )

## ----swiss_summary, eval=FALSE, echo=FALSE------------------------------------
# cat("=== SWISS ROLL MANIFOLD DATA ANALYSIS SUMMARY ===\n")
# cat("Dataset characteristics:\n")
# cat("- Objects:", swiss_results$data_characteristics$n_objects, "\n")
# cat("- Missing data:", round(swiss_results$data_characteristics$missing_percentage, 1), "%\n")
# cat("- Non-Euclidean deviation score:", round(swiss_results$data_characteristics$deviation_score, 4), "\n")
# cat("- Dimensions for 90% variance:", swiss_results$data_characteristics$dims_90_percent, "\n")
# cat("- Target embedding dimensions:", swiss_results$target_dims, "\n\n")
# 
# cat("Performance Summary:\n")
# if(sum(swiss_results$valid_topolow_results) > 0) {
#   cat("- Topolow mean error:", round(mean(swiss_results$topolow_errors[swiss_results$valid_topolow_results]), 2),
#       "±", round(sd(swiss_results$topolow_errors[swiss_results$valid_topolow_results]), 2), "\n")
# }
# if(!is.na(swiss_results$classical_mds_error)) {
#   cat("- Classical MDS error:", round(swiss_results$classical_mds_error, 2), "\n")
# }
# if(sum(swiss_results$valid_iterative_results) > 0) {
#   cat("- Iterative MDS mean error:", round(mean(swiss_results$iterative_mds_errors[swiss_results$valid_iterative_results]), 2),
#       "±", round(sd(swiss_results$iterative_mds_errors[swiss_results$valid_iterative_results]), 2), "\n")
# }

## ----swiss_eigenvalues, eval=FALSE, echo=FALSE--------------------------------
# # Create eigenvalue visualization for Swiss Roll data
# D_squared <- swiss_results$complete_distances^2
# n <- nrow(D_squared)
# J <- diag(n) - (1/n) * matrix(1, n, n)
# B <- -0.5 * J %*% D_squared %*% J
# eigenvals <- eigen(B, symmetric = TRUE)$values
# 
# eigenval_df <- data.frame(
#   index = 1:length(eigenvals),
#   eigenvalue = eigenvals,
#   type = ifelse(eigenvals > 1e-12, "Positive",
#                ifelse(eigenvals < -1e-12, "Negative", "Zero"))
# )
# 
# p_eigen_swiss <- ggplot(eigenval_df, aes(x = index, y = eigenvalue, fill = type)) +
#   geom_bar(stat = "identity", alpha = 0.8) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1, color = "black") +
#   scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "red", "Zero" = "gray")) +
#   labs(
#     title = "Eigenvalue Spectrum: Swiss Roll Manifold Data",
#     subtitle = paste("Deviation Score:", round(swiss_results$data_characteristics$deviation_score, 4),
#                      "| Non-Euclidean Variance:", round(swiss_results$data_characteristics$negative_variance_fraction * 100, 1), "%"),
#     x = "Eigenvalue Index (sorted descending)",
#     y = "Eigenvalue",
#     fill = "Eigenvalue Type"
#   ) +
#   theme_minimal() +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(size = 12, face = "bold"),
#     panel.grid.minor = element_blank()
#   )
# 
# print(p_eigen_swiss)

## ----swiss_eigenvalues-1, echo=FALSE, out.width="100%"------------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/swiss_eigenvalues-1.png")

## ----swiss_performance, eval=FALSE, echo=FALSE--------------------------------
# # Create performance comparison for Swiss Roll data
# results_df_swiss <- data.frame()
# 
# if(sum(swiss_results$valid_topolow_results) > 0) {
#   topolow_df <- data.frame(
#     Run = which(swiss_results$valid_topolow_results),
#     Method = "Topolow",
#     Total_Error = swiss_results$topolow_errors[swiss_results$valid_topolow_results],
#     Method_Type = "Stochastic",
#     Dataset = "Swiss Roll"
#   )
#   results_df_swiss <- rbind(results_df_swiss, topolow_df)
# }
# 
# if(!is.na(swiss_results$classical_mds_error)) {
#   classical_df <- data.frame(
#     Run = 1,
#     Method = "Classical MDS",
#     Total_Error = swiss_results$classical_mds_error,
#     Method_Type = "Deterministic",
#     Dataset = "Swiss Roll"
#   )
#   results_df_swiss <- rbind(results_df_swiss, classical_df)
# }
# 
# if(sum(swiss_results$valid_iterative_results) > 0) {
#   iterative_df <- data.frame(
#     Run = which(swiss_results$valid_iterative_results),
#     Method = "Iterative MDS",
#     Total_Error = swiss_results$iterative_mds_errors[swiss_results$valid_iterative_results],
#     Method_Type = "Stochastic",
#     Dataset = "Swiss Roll"
#   )
#   results_df_swiss <- rbind(results_df_swiss, iterative_df)
# }
# 
# if(nrow(results_df_swiss) > 0) {
#   p_performance_swiss <- ggplot(results_df_swiss, aes(x = Method, y = Total_Error, fill = Method_Type)) +
#     geom_boxplot(alpha = 0.7, outlier.shape = NA) +
#     geom_jitter(width = 0.2, alpha = 0.7, size = 2.5) +
#     scale_fill_manual(values = c("Stochastic" = "steelblue", "Deterministic" = "darkorange")) +
#     labs(
#       title = "Embedding Performance: Swiss Roll Manifold Data",
#       subtitle = paste("Non-Euclidean Score:", round(swiss_results$data_characteristics$deviation_score, 4),
#                        "| Missing Data:", round(swiss_results$data_characteristics$missing_percentage, 1), "%"),
#       y = "Sum of Absolute Distance Errors (L1 Norm)",
#       x = "Method",
#       fill = "Method Type"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 12, face = "bold"),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom"
#     )
# 
#   print(p_performance_swiss)
# }

## ----swiss_performance-1, echo=FALSE, out.width="100%"------------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/swiss_performance-1.png")

## ----combined_analysis, eval=FALSE, echo=FALSE--------------------------------
# # Combine results from both datasets
# combined_results <- rbind(
#   if(nrow(results_df_distorted) > 0) results_df_distorted else data.frame(),
#   if(nrow(results_df_swiss) > 0) results_df_swiss else data.frame()
# )
# 
# if(nrow(combined_results) > 0) {
#   # Overall performance comparison
#   p_combined_performance <- ggplot(combined_results, aes(x = Method, y = Total_Error)) +
#     geom_boxplot(alpha = 0.7, fill = "lightblue") +
#     geom_point(position = position_jitter(width = 0.2),
#                size = 2, alpha = 0.8, color = "darkblue") +
#     facet_wrap(~ Dataset, scales = "free_y") +
#     labs(
#       title = "Comparative Performance Across Dataset Types",
#       subtitle = "Error comparison for three methods on Synthesized non-Euclidean and Swiss Roll manifold data",
#       y = "Sum of Absolute Distance Errors (L1 Norm)",
#       x = "Method"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 12, face = "bold"),
#       panel.grid.minor = element_blank(),
#       strip.text = element_text(face = "bold", size = 10),
#       axis.text.x = element_text(angle = 45, hjust = 1)
#     )
# 
#   print(p_combined_performance)
# 
#   # Summary statistics by method and dataset
#   summary_stats <- combined_results %>%
#     group_by(Method, Dataset, Method_Type) %>%
#     summarise(
#       Count = n(),
#       Mean_Error = mean(Total_Error),
#       SD_Error = sd(Total_Error),
#       Min_Error = min(Total_Error),
#       Max_Error = max(Total_Error),
#       .groups = 'drop'
#     )
# 
#   cat("\n=== COMBINED PERFORMANCE SUMMARY ===\n")
#   print(summary_stats)
# }

## ----combined_analysis-1, echo=FALSE, out.width="100%"------------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/combined_analysis-1.png")

## ----statistical_analysis, eval=FALSE, echo=FALSE-----------------------------
# cat("\n=== STATISTICAL ANALYSIS ===\n")
# 
# # Statistical comparisons for Synthesized non-Euclidean data
# if(sum(distorted_results$valid_topolow_results) > 0 &&
#    sum(distorted_results$valid_iterative_results) > 0) {
#   cat("\nSynthesized non-Euclidean Data:\n")
#   enhanced_statistical_comparison(
#     distorted_results$topolow_errors[distorted_results$valid_topolow_results],
#     distorted_results$iterative_mds_errors[distorted_results$valid_iterative_results],
#     "Iterative MDS"
#   )
# }
# 
# # Statistical comparisons for Swiss Roll data
# if(sum(swiss_results$valid_topolow_results) > 0 &&
#    sum(swiss_results$valid_iterative_results) > 0) {
#   cat("\nSwiss Roll Manifold Data:\n")
#   enhanced_statistical_comparison(
#     swiss_results$topolow_errors[swiss_results$valid_topolow_results],
#     swiss_results$iterative_mds_errors[swiss_results$valid_iterative_results],
#     "Iterative MDS"
#   )
# }

## ----distance_preservation, eval=FALSE, echo=FALSE----------------------------
# # Distance preservation analysis for both datasets
# analyze_distance_preservation <- function(results, dataset_name) {
#   cat("\n=== DISTANCE PRESERVATION:", toupper(dataset_name), "===\n")
# 
#   # Get best performing runs
#   best_methods <- list()
# 
#   if(sum(results$valid_topolow_results) > 0) {
#     best_topolow_idx <- which.min(results$topolow_errors[results$valid_topolow_results])
#     actual_topolow_idx <- which(results$valid_topolow_results)[best_topolow_idx]
#     best_methods$topolow <- results$topolow_results[[actual_topolow_idx]]
#   }
# 
#   if(!is.na(results$classical_mds_error)) {
#     best_methods$classical <- results$classical_mds_result
#   }
# 
#   if(sum(results$valid_iterative_results) > 0) {
#     best_iterative_idx <- which.min(results$iterative_mds_errors[results$valid_iterative_results])
#     actual_iterative_idx <- which(results$valid_iterative_results)[best_iterative_idx]
#     best_methods$iterative <- results$iterative_mds_results[[actual_iterative_idx]]
#   }
# 
#   if(length(best_methods) > 0) {
#     # Calculate correlations
#     upper_tri_mask <- upper.tri(results$complete_distances)
#     true_distances_vector <- results$complete_distances[upper_tri_mask]
# 
#     distance_comparison_list <- list()
#     preservation_stats <- data.frame()
# 
#     for(method_name in names(best_methods)) {
#       embedded_distances <- as.matrix(dist(best_methods[[method_name]]))
#       embedded_distances_vector <- embedded_distances[upper_tri_mask]
# 
#       correlation <- cor(true_distances_vector, embedded_distances_vector, use = "complete.obs")
#       rsq <- summary(lm(embedded_distances_vector ~ true_distances_vector))$r.squared
# 
#       method_display_name <- switch(method_name,
#                                    "topolow" = "Topolow",
#                                    "classical" = "Classical MDS",
#                                    "iterative" = "Iterative MDS")
# 
#       preservation_stats <- rbind(preservation_stats, data.frame(
#         Method = method_display_name,
#         Dataset = dataset_name,
#         Correlation = correlation,
#         R_squared = rsq
#       ))
# 
#       distance_comparison_list[[method_name]] <- data.frame(
#         True_Distance = true_distances_vector,
#         Embedded_Distance = embedded_distances_vector,
#         Method = method_display_name,
#         Dataset = dataset_name
#       )
#     }
# 
#     distance_comparison_df <- do.call(rbind, distance_comparison_list)
# 
#     # Create distance preservation plot
#     p_distance_preservation <- ggplot(distance_comparison_df,
#                                      aes(x = True_Distance, y = Embedded_Distance, color = Method)) +
#       geom_point(alpha = 0.6, size = 1.5) +
#       geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
#       geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", alpha = 0.7) +
#       facet_wrap(~Method, scales = "free_x") +
#       coord_cartesian(ylim = c(0, NA)) +
#       labs(
#         title = paste("Distance Preservation:", dataset_name),
#         x = "Original Non-Euclidean Distance",
#         y = "Embedded Euclidean Distance"
#       ) +
#       theme_minimal() +
#       theme(
#         plot.title = element_text(size = 12, face = "bold"),
#         legend.position = "none",
#         panel.grid.minor = element_blank()
#       )
# 
#     print(p_distance_preservation)
# 
#     cat("\nDistance Preservation Statistics:\n")
#     print(preservation_stats)
# 
#     return(preservation_stats)
#   }
# }
# 
# # Analyze both datasets
# distorted_preservation <- analyze_distance_preservation(distorted_results, "Synthesized non-Euclidean")
# swiss_preservation <- analyze_distance_preservation(swiss_results, "Swiss Roll")

## ----distance_preservation-1, echo=FALSE, out.width="100%"--------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/distance_preservation-1.png")

## ----distance_preservation-2, echo=FALSE, out.width="100%"--------------------
# Include pre-generated figure from man/figures
knitr::include_graphics("../man/figures/distance_preservation-2.png")

## ----preservation_summary, eval=FALSE, echo=FALSE-----------------------------
# if(exists("distorted_preservation") && exists("swiss_preservation")) {
#   combined_preservation <- rbind(distorted_preservation, swiss_preservation)
# 
#   cat("Distance Preservation Summary (Correlation with Original Distances):\n")
#   print(combined_preservation)
# 
#   # Calculate mean correlations by method
#   mean_correlations <- combined_preservation %>%
#     group_by(Method) %>%
#     summarise(Mean_Correlation = mean(Correlation, na.rm = TRUE),
#               Mean_R_squared = mean(R_squared, na.rm = TRUE),
#               .groups = 'drop')
# 
#   cat("\nMean Performance Across Dataset Types:\n")
#   print(mean_correlations)
# }

## ----session_info, echo=FALSE-------------------------------------------------
sessionInfo()

