#' Multifile Record Linkage and Duplicate Detection
#'
#' The multilink package implements the methodology of Aleshin-Guendel & Sadinle
#' (2022). It handles the general problem of multifile record linkage and
#' duplicate detection, where any number of files are to be linked, and any of
#' the files may have duplicates.
#'
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\href{https://doi.org/10.1080/01621459.2021.2013242}{Published}] [\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#'
#' @examples
#' # Here we demonstrate an example workflow with the no duplicate dataset
#' data(no_dup_data)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(no_dup_data$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = no_dup_data$file_sizes,
#'  duplicates = c(0, 0, 0))
#'
#' # Specify the prior
#' prior_list <- specify_prior(comparison_list, mus = NA, nus = NA, flat = 0,
#'  alphas = rep(1, 7), dup_upper_bound = c(1, 1, 1),
#'  dup_count_prior_family = NA, dup_count_prior_pars = NA,
#'  n_prior_family = "uniform", n_prior_pars = NA)
#'
#' # Find initialization for the matching (this step is optional)
#' # The following line corresponds to only keeping pairs of records as
#' # potential matches in the initialization for which neither gname nor fname
#' # disagree at the highest level
#' pairs_to_keep <- (comparison_list$comparisons[, "gname_DL_3"] != TRUE) &
#'  (comparison_list$comparisons[, "fname_DL_3"] != TRUE)
#' Z_init <- initialize_partition(comparison_list, pairs_to_keep, seed = 42)
#'
#' # Run the Gibbs sampler
#' # Takes around 2 and a half minutes
#' \dontrun{
#' results <- gibbs_sampler(comparison_list, prior_list, n_iter = 1000,
#'  Z_init = Z_init, seed = 42)
#'
#' # Find the full Bayes estimate
#' # Takes a couple of seconds
#' full_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = Inf, max_cc_size = 50)
#'
#' # There are 492 clusters in the full estimate
#' length(unique(full_estimate))
#' # There are 500 entities represented in the records
#' length(unique(no_dup_data$IDs))
#'
#' # Find which record pairs are truly coreferent based on IDs
#' true_links <- no_dup_data$IDs[comparison_list$record_pairs[, 1]] ==
#' no_dup_data$IDs[comparison_list$record_pairs[, 2]]
#'
#' # Find which record pairs are in the same clusters in the full estimate
#' full_estimate_links <- full_estimate[comparison_list$record_pairs[, 1]] ==
#' full_estimate[comparison_list$record_pairs[, 2]]
#'
#' # Find the number of true matches in the full estimate
#' true_matches <- sum(full_estimate_links & true_links)
#'
#' # Precision of the full estimate is 0.96
#' true_matches / sum(full_estimate_links)
#'
#' # Recall of the full estimate is 1
#' true_matches / sum(true_links)
#'
#' # Find the partial Bayes estimate
#' # Takes around 20 seconds
#' partial_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = 0.1, max_cc_size = 12)
#'
#' # The partial estimate abstains from making decisions for 28 records
#' sum(partial_estimate == -1)
#'
#' # For the records which decisions were made for in the partial estimate,
#' # there are 473 clusters
#' length(unique(partial_estimate))
#'
#' # Abstain rate of partial_estimate is 0.04
#' sum(partial_estimate == -1) / length(partial_estimate)
#'
#' # Relabel records where we abstained
#' partial_estimate[which(partial_estimate == -1)] <- length(partial_estimate) +
#' which(partial_estimate == -1)
#'
#' # Find which record pairs are in the same clusters in the full estimate
#' partial_estimate_links <-
#'  partial_estimate[comparison_list$record_pairs[, 1]] ==
#'  partial_estimate[comparison_list$record_pairs[, 2]]
#'
#' # Find the number of true matches in the partial estimate
#' true_matches_A <- sum(partial_estimate_links & true_links)
#'
#' # Precision of the partial estimate is 0.98
#' true_matches_A / sum(partial_estimate_links)
#' }
#'
#' # Here we demonstrate an example workflow with the duplicate dataset
#' data(dup_data)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(dup_data$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = dup_data$file_sizes,
#'  duplicates = c(1, 1, 1))
#'
#' # Reduce the comparison data
#' # The following line corresponds to only keeping pairs of records for which
#' # neither gname nor fname disagree at the highest level
#' pairs_to_keep <- (comparison_list$comparisons[, "gname_DL_3"] != TRUE) &
#'  (comparison_list$comparisons[, "fname_DL_3"] != TRUE)
#' reduced_comparison_list <- reduce_comparison_data(comparison_list,
#'  pairs_to_keep, cc = 1)
#'
#' # Specify the prior
#' prior_list <- specify_prior(reduced_comparison_list, mus = NA, nus = NA,
#'  flat = 0, alphas = rep(1, 7), dup_upper_bound = c(10, 10, 10),
#'  dup_count_prior_family = c("Poisson", "Poisson", "Poisson"),
#'  dup_count_prior_pars = list(c(1), c(1), c(1)), n_prior_family = "uniform",
#'  n_prior_pars = NA)
#'
#' # Run the Gibbs sampler
#' # Takes around 10 seconds
#' \dontrun{
#' results <- gibbs_sampler(reduced_comparison_list, prior_list, n_iter = 1000,
#'  seed = 42)
#'
#' # Find the full Bayes estimate
#' # Takes a couple of seconds
#' full_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = Inf, max_cc_size = 50)
#'
#' # There are 495 clusters in the full estimate (including records records
#' # determined not to be candidate matches to any other records using
#' # reduce_comparison_data)
#' length(unique(full_estimate)) +
#' sum(reduced_comparison_list$file_sizes_not_included)
#' # There are 500 entities represented in the records
#' length(unique(dup_data$IDs))
#'
#' # Find which record pairs are truly coreferent based on IDs
#' true_links <- dup_data$IDs[comparison_list$record_pairs[, 1]] ==
#' dup_data$IDs[comparison_list$record_pairs[, 2]]
#'
#' # Focus on the record pairs that were candidate matches
#' true_links_reduced <- true_links[reduced_comparison_list$pairs_to_keep]
#'
#' # Calculate the number of prior false non-matches based on the indexing
#' # scheme used
#' prior_fnm <-
#'  nrow(comparison_list$record_pairs[true_links &
#'  (!reduced_comparison_list$pairs_to_keep), ])
#'
#' # Find which record pairs are in the same clusters in the full estimate
#' full_estimate_links <-
#'  full_estimate[reduced_comparison_list$record_pairs[, 1]] ==
#'  full_estimate[reduced_comparison_list$record_pairs[, 2]]
#'
#' # Find the number of true matches in the full estimate
#' true_matches <- sum(full_estimate_links & true_links_reduced)
#'
#' # Precision of the full estimate is 0.95
#' true_matches / sum(full_estimate_links)
#'
#' # Recall of the full estimate is 0.99
#' true_matches / (sum(true_links_reduced) + prior_fnm)
#'
#' # Find the partial Bayes estimate
#' # Takes around 15 seconds
#' partial_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = 0.1, max_cc_size = 12)
#'
#' # The partial estimate abstains from making decisions for 3 records
#' sum(partial_estimate == -1)
#'
#' # For the records which decisions were made for in the partial estimate,
#' # there are 495 clusters (including records determined not to be candidate
#' # matches to any other records using reduce_comparison_data)
#' length(unique(partial_estimate)) +
#'  sum(reduced_comparison_list$file_sizes_not_included)
#'
#' # Abstain rate of partial_estimate is 0.004 (excluding records determined not
#' # to be candidate matches to any other records using reduce_comparison_data)
#' sum(partial_estimate == -1) / length(partial_estimate)
#'
#' # Relabel records where we abstained
#' partial_estimate[which(partial_estimate == -1)] <- length(partial_estimate) +
#' which(partial_estimate == -1)
#'
#' # Find which record pairs are in the same clusters in the full estimate
#' partial_estimate_links <-
#'  partial_estimate[reduced_comparison_list$record_pairs[, 1]] ==
#'  partial_estimate[reduced_comparison_list$record_pairs[, 2]]
#'
#' # Find the number of true matches in the partial estimate
#' true_matches_A <- sum(partial_estimate_links & true_links_reduced)
#'
#' # Precision of the partial estimate is 0.96
#' true_matches_A / sum(partial_estimate_links)
#'
#' # Relabel the full and partial Bayes estimates
#' full_estimate_relabel <- relabel_bayes_estimate(reduced_comparison_list,
#'  full_estimate)
#'
#'  partial_estimate_relabel <- relabel_bayes_estimate(reduced_comparison_list,
#'   partial_estimate)
#'
#'  # Add columns to the records corresponding to their full and partial
#'  # Bayes estimates
#'  dup_data$records <- cbind(dup_data$records,
#'   full_estimate_id = full_estimate_relabel$link_id,
#'   partial_estimate_id = partial_estimate_relabel$link_id)
#' }
#'
#'
#' @docType package
#' @name multilink
#' @useDynLib multilink, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#> NULL
