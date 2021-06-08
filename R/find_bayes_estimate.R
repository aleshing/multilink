#' Find the Bayes Estimate of a Partition
#'
#' Find the (approximate) Bayes estimate of a partition based on MCMC samples
#' of the partition and a specified loss function.
#'
#' @param partitions Posterior samples of the partition, where each column
#' is one sample and the partition is represented as an \code{integer} vector of
#' arbitrary labels, as produced by the output of a call to
#' \code{\link{gibbs_sampler}}.
#' @param burn_in The number of samples to discard for burn in.
#' @param L_FNM Positive loss for a false non-match. Default is \code{1}.
#' @param L_FM1 Positive loss for a type 1 false match. Default is \code{1}.
#' @param L_FM2 Positive loss for a type 2 false match. Default is \code{2}.
#' @param L_A Positive loss for abstaining from making a decision for a record.
#' Default is \code{Inf}, i.e. decisions are made for all records.
#' @param max_cc_size The maximum allowable connected component size over which
#' the posterior expected loss is minimized. Default is \code{nrow(partitions)},
#' i.e. no approximation is used. When \code{is.infinite(L_A)}, we recommend
#' setting this argument to \code{50}, then increasing based on a computational
#' budget (setting this argument to values around \code{100-200} should be
#' computationally feasible). When \code{!is.infinite(L_A)}, we recommend
#' setting this argument to \code{10-12}, then increasing based on a
#' computational budget (although an increase of \code{1} in this argument can
#' in the worst case lead to a doubling in computation time).
#'
#' @return
#' A vector, the same length of a column of \code{partitions} containing the
#' (approximate) Bayes estimate of the partition. If \code{!is.infinite(L_A)}
#' the output may be a partial estimate. A positive number \code{l} in index
#' \code{i} indicates that record \code{i} is in the same cluster as every other
#' record \code{j} with \code{l} in index \code{j}. A value of \code{-1} in
#' index \code{i} indicates that the Bayes estimate abstained from  making a
#' decision for record \code{i}.
#' @export
#'
#' @examples
#' # Example with no duplicate dataset
#' data(no_dup_data)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(no_dup_data$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5),  c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = no_dup_data$file_sizes,
#'  duplicates = c(0, 0, 0))
#'
#' # Specify the prior
#' prior_list <- specify_prior(comparison_list, mus = NA, nus = NA, flat = 0,
#'  alphas = rep(1, 7), dup_upper_bound = c(1, 1, 1),
#'  dup_count_prior_family = NA, dup_count_prior_pars = NA,
#'  n_prior_family = "uniform", n_prior_pars = NA)
#'
#' # Find initialization for the matching
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
#' # Find the partial Bayes estimate
#' # Takes around 20 seconds
#' partial_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = 0.1, max_cc_size = 12)
#' }
#'
#' # Example with duplicate dataset
#' data(dup_data)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(dup_data$records,
#' types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
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
#'  dup_count_prior_pars = list(c(1), c(1), c(1)),
#'  n_prior_family = "uniform", n_prior_pars = NA)
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
#' # Find the partial Bayes estimate
#' # Takes around 15 seconds
#' partial_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = 0.1, max_cc_size = 12)
#' }
find_bayes_estimate <- function(partitions, burn_in, L_FNM = 1, L_FM1 = 1,
                                L_FM2 = 2, L_A = Inf,
                                max_cc_size = nrow(partitions)){
    # Input checks
    if(!is.numeric(L_FNM)){
        stop("'L_FNM' must be numeric")
    }
    if(L_FNM <= 0){
        stop("'L_FNM' must be positive")
    }
    if(!is.numeric(L_FM1)){
        stop("'L_FM1' must be numeric")
    }
    if(L_FM1 <= 0){
        stop("'L_FM1' must be positive")
    }
    if(!is.numeric(L_FM2)){
        stop("'L_FM2' must be numeric")
    }
    if(L_FM2 <= 0){
        stop("'L_FM2' must be positive")
    }
    if(!is.numeric(L_A)){
        stop("'L_A' must be numeric")
    }
    if(L_A <= 0){
        stop("'L_A' must be positive")
    }
    if(!is.infinite(L_A) & max_cc_size > 15){
        print("Finding partial Bayes estimate may take a long time when maximum
              connected component is size > 15")
    }
    if(is.infinite(L_A) & max_cc_size > 100){
        print("Finding full Bayes estimate may take a long time when maximum
              connected component is size > 100")
    }

    n_iter <- ncol(partitions)
    partitions <- as.matrix(partitions[, (burn_in + 1):n_iter])
    r <- nrow(partitions)
    TT <- n_iter - burn_in

    thresh <- 0
    thresh_inc <- 1 / TT

    temp_psm <- mcclust::comp.psm(t(partitions))
    psm <- temp_psm > thresh
    g <- igraph::graph_from_adjacency_matrix(psm, mode = "undirected")
    ccs <- igraph::components(g)

    while(max(ccs$csize) > max_cc_size){
        thresh <- thresh + thresh_inc
        psm <- temp_psm > thresh
        g <- igraph::graph_from_adjacency_matrix(psm, mode = "undirected")
        ccs <- igraph::components(g)
    }

    print(paste0("Finding Bayes estimate with a threshold of ", thresh,
                 " and a maximum connected component of size ", max(ccs$csize)))

    relabel_partition <- function(sample){
        iterlabels <- unique(sample)
        newlabs <- seq_len(length(iterlabels))
        newlabs[match(sample, iterlabels)]
    }

    Z_hat <- rep(0, r)
    for(cc in 1:ccs$no){
        temp <- which(ccs$membership == cc)
        if(ccs$csize[cc] == 1){
            Z_hat[temp] <- max(Z_hat) + 1
        }
        else{
            temp_part <- partitions[temp, ]
            if(is.infinite(L_A)){
                temp_loss <- get_posterior_loss_allcpp(TT, nrow(temp_part),
                                                       temp_part, L_FNM, L_FM1,
                                                       L_FM2)
                Z_hat[temp] <-
                    relabel_partition(temp_part[, which.min(temp_loss)]) +
                    max(Z_hat)
            }
            else{
                power <- as.matrix(expand.grid(replicate(ccs$csize[cc], c(1, 0),
                                                         simplify = FALSE)))
                temp_loss <- get_posterior_loss_abstain_cpp(TT, nrow(temp_part),
                                                            temp_part, L_FNM,
                                                            L_FM1, L_FM2, L_A,
                                                            power)
                min_loss <- which.min(temp_loss)
                min_loss_draw <- ((min_loss - 1) %% TT) + 1
                min_loss_power_row <- ((min_loss - 1) %/% TT) + 1
                temp_Z <- temp_part[, min_loss_draw]
                temp_decision <- which(power[min_loss_power_row, ] == 1)
                temp_Z[temp_decision] <-
                    relabel_partition(temp_Z[temp_decision]) + max(Z_hat)
                temp_Z[which(power[min_loss_power_row,] == 0)] <- -1
                Z_hat[temp] <- temp_Z
            }
        }
    }

    return(Z_hat)
}
