#' Gibbs Sampler for Posterior Inference
#'
#' Run a Gibbs sampler to explore the posterior distribution of partitions of
#' records.
#'
#' Given the prior specified using \code{\link{specify_prior}}, this function
#' runs a Gibbs sampler to explore the posterior distribution of partitions of
#' records, conditional on the comparison data created using
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}.
#'
#' @param comparison_list The output from a call to
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}.
#' @param prior_list The output from a call to \code{\link{specify_prior}}.
#' @param n_iter The number of iterations of the Gibbs sampler to run.
#' @param Z_init Initialization of the partition of records, represented as an
#' \code{integer} vector of arbitrary labels of length
#' \code{sum(comparison_list$file_sizes)}. The default initialization places
#' each record in its own cluster. See \code{\link{initialize_partition}} for an
#' alternative initialization when there are no duplicates in each file.
#' @param seed The seed to use while running the Gibbs sampler.
#' @param single_likelihood A \code{logical} indicator of whether to use a
#' single likelihood for comparisons for all file pairs, or whether to use a
#' separate likelihood for comparisons for each file pair. When
#' \code{single_likelihood=TRUE}, a single likelihood is used, and the prior
#' hyperparameters for \code{m} and \code{u} from the first file pair are used.
#' We do not recommend using a single likelihood in general.
#'
#' @return a list containing:
#' \describe{
#'   \item{\code{m}}{Posterior samples of the \code{m} parameters. Each column
#'   is one sample.}
#'   \item{\code{u}}{Posterior samples of the \code{u} parameters. Each column
#'   is one sample.}
#'   \item{\code{partitions}}{Posterior samples of the partition. Each column
#'   is one sample. Note that the partition is represented as an \code{integer}
#'   vector of arbitrary labels of length
#'   \code{sum(comparison_list$file_sizes)}.}
#'   \item{\code{contingency_tables}}{Posterior samples of the overlap table.
#'   Each column is one sample. This incorporates counts of records determined
#'   not to be candidate matches to any other records using
#'   \code{\link{reduce_comparison_data}}.}
#'   \item{\code{cluster_sizes}}{Posterior samples of the size of each cluster
#'   (associated with an arbitrary label from \code{1} to
#'   \code{sum(comparison_list$file_sizes)}). Each column is one sample.}
#'   \item{\code{sampling_time}}{The time in seconds it took to run the
#'   sampler.}
#' }
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
#' }
gibbs_sampler <- function(comparison_list, prior_list, n_iter = 2000,
                          Z_init = 1:sum(comparison_list$file_sizes),
                          seed = 70, single_likelihood = FALSE){
    # Input checks
    if(length(Z_init) != sum(comparison_list$file_sizes)){
        stop("'length(Z_init)' and 'sum(comparison_list$file_sizes)' must be the
             same")
    }
    if(sum(comparison_list$ab_not_included) > 0 &
       sum(prior_list$nus - prior_list$nus_specified > 0) == 0){
        stop("Check to see if comparison_list is the same as comparison_list
                used in specify_prior.")
    }
    if(sum(comparison_list$ab_not_included) == 0 &
       sum(prior_list$nus - prior_list$nus_specified > 0) > 0){
        stop("Check to see if comparison_list is the same as comparison_list
                used in specify_prior.")
    }

    # Progress
    print("Processing inputs")

    # Dump the contents of comparison_list
    # Converted record_pairs to a matrix for Rcpp
    record_pairs <- as.matrix(comparison_list$record_pairs)
    comparisons <- comparison_list$comparisons
    K <- comparison_list$K
    file_sizes <- comparison_list$file_sizes
    duplicates <- comparison_list$duplicates
    field_levels <- comparison_list$field_levels
    file_labels <- comparison_list$file_labels
    fp_matrix <- comparison_list$fp_matrix
    rp_to_fp <- comparison_list$rp_to_fp
    ab <- comparison_list$ab
    file_sizes_not_included <- comparison_list$file_sizes_not_included
    ab_not_included <- comparison_list$ab_not_included
    cc <- comparison_list$cc


    # A hack to overcome how I initialize the partition without initializing the
    # model parameters
    n_iter <- n_iter + 1

    # r: the total number of records
    r <- sum(file_sizes)

    # Z_samp: a matrix of size r x n_iter, where Z_samp[, i] is the partition at
    # iteration i, using an arbitrary labeling (where the labels take values in
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

    # powers: Given a binary vector representing the index of a cell in a 2 ^ K
    # contingency table, powers is used to to convert the index to an integer
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
    # labeling of the rows is the same labeling as for cont
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

    # Dump the contents of prior_list
    flat <- prior_list$flat
    mus <- prior_list$mus
    nus <- prior_list$nus
    nus_specified <- prior_list$nus_specified
    no_dups  <- prior_list$no_dups
    alphas <- prior_list$alphas
    alpha_0 <- prior_list$alpha_0
    dup_upper_bound <- prior_list$dup_upper_bound
    log_dup_count_prior <- prior_list$log_dup_count_prior
    log_n_prior <- prior_list$log_n_prior

    # Determine whether indexing was used
    indexing_used <- 0
    if(sum(ab_not_included) > 0){
        indexing_used <- 1
    }

    # Calculate relevant summaries if only a single likelihood is used
    single_nus <- rep(0, L)
    single_ab <- rep(0, L)
    if(single_likelihood){
        single_nus <- nus_specified[1:L]
        for(f in 1:num_fp){
            single_nus <- single_nus +
                ab_not_included[((f - 1) * L + 1):(f * L)]
            single_ab <- single_ab + ab[((f - 1) * L + 1):(f * L)]
        }
    }

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
                                    cc, Z_members, clust_sizes_collapsed,
                                    indexing_used, single_likelihood,
                                    single_nus, single_ab)

    # Keep track of time
    end <- Sys.time()
    sampling_time <- as.numeric(end-start, units = "secs")
    # Progress
    print(paste0("Finished sampling in ", sampling_time, " seconds"))

    # Return the posterior samples and the derived contingency tables and
    # cluster sizes, removing the first column
    contingency_tables <- gibbs_output$cont_samp[, 2:n_iter]
    if(!single_likelihood){
        m <- gibbs_output$m_samp[, 2:n_iter]
        u <- gibbs_output$u_samp[, 2:n_iter]
        rownames(m) <- names(ab)
        rownames(u) <- names(ab)
    }
    else{
        m <- gibbs_output$m_samp[1:ncol(comparisons), 2:n_iter]
        u <- gibbs_output$u_samp[1:ncol(comparisons), 2:n_iter]
        rownames(m) <- colnames(comparisons)
        rownames(u) <- colnames(comparisons)
    }
    rownames(contingency_tables) <- names(alphas)
    return(list(m = m, u = u, partitions = gibbs_output$Z_samp[, 2:n_iter],
                contingency_tables = contingency_tables,
                cluster_sizes = gibbs_output$clust_sizes_samp[, 2:n_iter],
                sampling_time = sampling_time))

}
