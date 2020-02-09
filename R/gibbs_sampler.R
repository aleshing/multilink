#' Gibbs Sampler for Posterior Inference
#'
#' Description goes here.
#'
#' @param comparison_list the output from a call to
#' \code{create_comparison_data} or \code{reduce_comparison_data}
#' @param partition_prior_list the output of a call to
#' \code{specify_partition_prior}
#' @param n_iter The number of iterations of the Gibbs sampler to run
#' @param mus,nus
#' @param Z_init
#' @param seed The seed to use while running the Gibbs sampler
#'
#' @return a list containing:
#' \describe{
#'   \item{\code{match_samp}}{}
#'   \item{\code{non_match_samp}}{}
#'   \item{\code{partitions}}{}
#'   \item{\code{contingency_tables}}{}
#'   \item{\code{cluster_sizes}}{}
#'   \item{\code{sampling_time}}{}
#' }
#' @export
#'
#' @examples
gibbs_sampler <- function(comparison_list, partition_prior_list, n_iter = 1000,
                          mus = NA, nus = NA, Z_init = 1:r, seed = 70){

    # Progress
    print("Processing inputs")

    # Dump the contents of comparison_list
    # Converted to a matrix for Rcpp
    record_pairs <- as.matrix(comparison_list$record_pairs)
    comparisons <- comparison_list$comparisons
    K <- comparison_list$K
    file_sizes <- comparison_list$file_sizes
    duplicates <- comparison_list$duplicates
    field_levels <- comparison_list$field_levels
    file_labels <- comparison_list$file_labels
    rp_to_fp <- comparison_list$rp_to_fp
    ab <- comparison_list$ab
    file_sizes_not_included <- comparison_list$file_sizes_not_included
    ab_not_included <- comparison_list$ab_not_included
    cc <- comparison_list$cc

    # fp_matrix <- create_fp_matrix(K)
    temp_fp <- expand.grid(1:K, 1:K)[2:1]
    temp_fp <- temp_fp[temp_fp[, 1] <= temp_fp[, 2], ]
    fp_matrix <- matrix(0, nrow = K, ncol = K)
    for(i in 1:nrow(temp_fp)){
        fp_matrix[temp_fp[i, 1], temp_fp[i, 2]] <- i
        fp_matrix[temp_fp[i, 2], temp_fp[i, 1]] <- i
    }

    # A hack to overcome how I initialize the partition without initializing the
    # model parameters
    n_iter <- n_iter + 1

    # r: the total number of records
    r <- sum(file_sizes)

    # Z_samp: a matrix of size r x n_iter, where Z_samp[, i] is the partition at
    # iteration i, using an arbitrary labelling (where the labels take values in
    # 1:r)
    Z_samp <- matrix(NA, nrow = r, ncol = n_iter)
    Z_samp[, 1] <- Z_init

    # Z_members: a matrix of size r x r, where the row Z_mat[i, ] contains the
    # records assigned to label i, followed by 0s
    Z_members <- matrix(0, nrow = r, ncol = r)
    if(sum(Z_init == 1:r) == r){ Z_members[, 1] <- 1:r }
    else{
        for(i in 1:max(Z_init)){
            members <- which(Z_init == i)
            Z_members[i, 1:length(members)] <- members
        }
    }

    # Generate powers for converting binary vectors to integers
    powers <- rev(2 ^ (0:(K - 1)))

    # cont: a vector of length 2^K representing the overlap table derived from
    # the current partition. Note that cont[i] corresponds to the cell
    # count of the overlap table for the inclusion pattern given by the binary
    # representation of i-1 with K bits
    cont <- rep(0, 2 ^ K)

    # clust_sizes: a matrix of size K x r, where clust_sizes[k, i] represents
    # the number of records from file k in the cluster with label i in the
    # current partition
    # n: the number of entities derived from the current partition
    clust_sizes <- matrix(0, nrow = K, ncol = r)
    if(sum(Z_init == 1:r) == r){
        for(i in 1:r){ clust_sizes[file_labels[i], i] <- 1 }
        # Initialize the overlap table to the partition of singleton clusters
        for(i in 1:r){
            ind <- sum(as.integer(clust_sizes[, i] > 0) * powers) + 1
            cont[ind] <- cont[ind] + 1
        }
        n <- r
    }
    else{
        for(i in 1:max(Z_init)){
            members <- which(Z_init == i)
            files <- file_labels[members]
            for(k in 1:K){ clust_sizes[k, i] <- sum(files == k) }
        }
        # Initialize the overlap table to that induced by Z_init
        for(i in 1:r){
            ind <- sum(as.integer(clust_sizes[, i] > 0) * powers) + 1
            if(ind > 1){ cont[ind] <- cont[ind] + 1 }
        }
        n <- length(unique(Z_init))
    }


    # clust_sizes_collapsed: a vector of size r, where clust_sizes[i] represents
    # the number of records in the cluster with label i in the current partition
    clust_sizes_collapsed <- colSums(clust_sizes)

    # singleton_ind: a length K vector, where singleton_ind[k] is the index in
    # the overlap table for a cluster with only records from file k
    singleton_ind <- rep(0, K)
    for(k in 1:K){
        ind <- rep(0, K)
        ind[k] <- 1
        singleton_ind[k] <- sum(ind * powers) + 1
    }

    # Incorporate the records not included in the linkage into the overlap table
    # and n.
    for(k in 1:K){
        cont[singleton_ind[k]] <- cont[singleton_ind[k]] +
            file_sizes_not_included[k]
    }
    n <- n + sum(file_sizes_not_included)

    # cont_samp: a matrix of size 2^K x n_iter, where cont_samp[, i] is the
    # overlap table derived from the partition at iteration i. Note that the
    # labelling of the rows is the same labelling as for cont
    cont_samp <- matrix(NA, nrow = 2^K, ncol = n_iter)
    cont_samp[, 1] <- cont

    # clust_sizes_samp: a matrix of size r x n_iter, where
    # clust_sizes_samp[i, j] represents the number of records in the cluster
    # with label j in iteration i
    clust_sizes_samp <- matrix(NA, nrow = r, ncol = n_iter)
    clust_sizes_samp[, 1] <- clust_sizes_collapsed

    # num_fp: the number of file pairs
    num_fp <- K * (K + 1)/2

    # FF: the number of fields being compared
    FF <- length(field_levels)

    # L: the sum of the number of agreement categories for each field
    L <- sum(field_levels)

    # level_cum: level_cum[f] gives the start of field f's levels
    # for easy indexing of m_samp and u_samp
    level_cum <- c(0, cumsum(field_levels)[1:FF]) + 1

    # num_rp: the number of record pairs we have comparison data for
    num_rp <- nrow(record_pairs)

    # valid_fp: a vector containing the integer indexes for file pairs being
    # compared, in C++ indexing
    temp_valid_fp <- upper.tri(fp_matrix)
    for(k in 1:K){ temp_valid_fp[k, k] <- as.logical(duplicates[k])}
    valid_fp <- fp_matrix[temp_valid_fp] - 1

    # rp_ind: the indices of record pairs in the comparison data. rp_ind[i, j]
    # gives the row index of comp_data where records i and j are compared, or -1
    # if the records aren't being compared
    rp_ind <- matrix(-1, nrow = r, ncol = r)
    for(rp in 1:num_rp){
        i <- record_pairs[rp, 1]
        j <- record_pairs[rp, 2]
        rp_ind[i, j] <- rp
        rp_ind[j, i] <- rp
    }

    # valid_rp: indicators of which pairs of records are possible links.
    # rp_ind[i, j] is TRUE if records i and j are compared (or if i=j, which
    # makes finding valid clusters easier), and FALSE if the records aren't
    # being compared
    valid_rp <- rp_ind > 0
    for(i in 1:r){ valid_rp[i, i] <- TRUE}

    # If no values are given for mus or nus, use default flat prior
    if(is.na(mus)){ mus <- rep(1, num_fp * L) }
    if(is.na(nus)){ nus <- rep(1, num_fp * L) }

    #### Add in prior non-coreferent counts using ab_full
    nus <- nus + ab_not_included

    # m_samp: the samples for the match model probabilities, where each column
    # is the sample from an iteration, structured the same way as mus/nus
    # u_samp: defined similarly to m_samp, but for non-match model probabilities
    m_samp <- matrix(0, nrow = num_fp * L, ncol = n_iter)
    u_samp <- matrix(0, nrow = num_fp * L, ncol = n_iter)

    # NOTE: m_samp/u_samp don't need to be initialized if Z_samp is initialized
    # and m_samp/u_samp are sampled first!

    # r_1: The index of the first record to sample a cluster membership for. In
    # particular, this is just 0 if you allow duplicates. If you assume any of
    # your files contains no duplicates, and you encode this in the
    # dup_upper_bound input as 1, we assume the largest such file with no
    # duplicates is file 1, and you won't sample cluster memberships for records
    # in file 1.
    r_1 <- 0
    if(duplicates[1] == 0){ r_1 <- file_sizes[1] }

    # Dump the contents of partition_prior_list
    flat <- partition_prior_list$flat
    no_dups  <- partition_prior_list$no_dups
    alphas <- partition_prior_list$alphas
    alpha_0 <- partition_prior_list$alpha_0
    dup_upper_bound <- partition_prior_list$dup_upper_bound
    log_dup_count_prior <- partition_prior_list$log_dup_count_prior
    log_n_prior <- partition_prior_list$log_n_prior

    set.seed(seed)

    # Progress
    print("Beginning sampling")

    # Keep track of time
    start <- Sys.time()

    gibbs_output <- gibbs_loop_rcpp(n_iter, Z_samp, clust_sizes_samp, cont_samp,
                                    m_samp, u_samp, mus, nus, alphas, alpha_0,
                                    dup_upper_bound, log_dup_count_prior,
                                    log_n_prior, cont, clust_sizes, n, ab,
                                    comparisons, record_pairs, flat, r, r_1,
                                    valid_rp, singleton_ind, rp_ind,
                                    file_labels, powers, L, num_fp, num_rp,
                                    FF, rp_to_fp, level_cum, no_dups, valid_fp,
                                    cc, Z_members, clust_sizes_collapsed)

    # Keep track of time
    end <- Sys.time()
    sampling_time <- as.numeric(end-start, units="secs")

    # Progress
    print(paste0("Finished sampling in ", sampling_time, " seconds"))

    # Return the posterior samples and the derived contingency tables and
    # cluster sizes, removing the first column
    return(list(match_samp=gibbs_output$m_samp[, 2:n_iter],
                non_match_samp=gibbs_output$u_samp[, 2:n_iter],
                partitions=gibbs_output$Z_samp[, 2:n_iter],
                contingency_tables=gibbs_output$cont_samp[, 2:n_iter],
                cluster_sizes=gibbs_output$clust_sizes_samp[, 2:n_iter],
                sampling_time=sampling_time))

}
