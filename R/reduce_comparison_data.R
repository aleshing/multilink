#' Reduce Comparison Data Size
#'
#' Description goes here.
#'
#' @param comparison_list the output from a call to
#' \code{\link{create_comparison_data}}
#' @param pairs_to_keep a logical vector indicating which record pairs should be
#' kept
#' @param cc an indicator of whether we should use the connected components
#' induced by pairs_to_keep
#'
#' @return a list containing:
#' \describe{
#'   \item{\code{record_pairs}}{}
#'   \item{\code{comparisons}}{}
#'   \item{\code{K}}{}
#'   \item{\code{file_sizes}}{}
#'   \item{\code{duplicates}}{}
#'   \item{\code{field_levels}}{}
#'   \item{\code{file_labels}}{}
#'   \item{\code{rp_to_fp}}{}
#'   \item{\code{ab}}{}
#'   \item{\code{file_sizes_not_included}}{}
#'   \item{\code{ab_not_included}}{}
#'   \item{\code{labels}}{}
#'   \item{\code{pairs_to_keep}}{}
#'   \item{\code{cc}}{}
#' }
#' @export
#'
#' @examples
#'
reduce_comparison_data <- function(comparison_list, pairs_to_keep, cc = 1){
    # Dump the contents of comparison_list
    # record_pairs <- comparison_list$record_pairs
    # comparisons <- comparison_list$comparisons
    # K <- comparison_list$K
    # file_sizes <- comparison_list$file_sizes
    # duplicates <- comparison_list$duplicates
    # field_levels <- comparison_list$field_levels
    # file_labels <- comparison_list$file_labels
    # fp_matrix <- create_fp_matrix(K)
    # rp_to_fp <- comparison_list$rp_to_fp
    # ab <- comparison_list$ab
    # file_sizes_not_included <- comparison_list$file_sizes_not_included
    # ab_not_included <- comparison_list$ab_not_included

    record_pairs <- comparison_list$record_pairs
    comparisons <- comparison_list$comparisons
    file_sizes <- comparison_list$file_sizes
    K <- comparison_list$K
    field_levels <- comparison_list$field_levels

    if(cc == 1){
        g <- igraph::graph_from_edgelist(
            as.matrix(record_pairs[pairs_to_keep, ]), directed=FALSE)
        components <- igraph::components(g)

        new_kept_pairs <- data.frame("i" = NULL, "j" = NULL)
        for(i in 1:length(components$csize)){
            if(components$csize[i] > 1){
                ind <- which(components$membership == i)
                # Note we don't need to construct these connected components
                # carefully (the data.table code deals with this)
                new_kept_pairs <- rbind(new_kept_pairs, t(utils::combn(ind, 2)))
            }
        }

        colnames(new_kept_pairs) <- c("i", "j")
        colnames(record_pairs) <- c("i", "j")

        # record_pairs_dt <- data.table::data.table(
        #     cbind(record_pairs, "record_pairs_rows" = 1:nrow(record_pairs)))
        # data.table::setkey(record_pairs_dt, 'i', 'j')
        # # record_pairs_dt[, "record_pairs_rows" := .I]
        # new_kept_pairs_dt <- data.table::data.table(new_kept_pairs)
        # data.table::setkey(new_kept_pairs_dt, 'i', 'j')
        # transitive_closure <- record_pairs_dt[new_kept_pairs_dt]
        # new_pairs_to_keep <- rep(FALSE, nrow(record_pairs))
        # new_pairs_to_keep[transitive_closure$record_pairs_rows] <- TRUE

        record_pairs_ids <- paste(record_pairs[, 1], record_pairs[, 2])
        new_kept_pairs_ids <- paste(new_kept_pairs[, 1], new_kept_pairs[, 2])
        common_ids <- which(record_pairs_ids %in% new_kept_pairs_ids)
        new_pairs_to_keep <- rep(FALSE, nrow(record_pairs))
        new_pairs_to_keep[common_ids] <- TRUE
        pairs_to_keep <- new_pairs_to_keep
    }

    kept_pairs <- record_pairs[pairs_to_keep, ]
    kept_comparisons <- comparisons[pairs_to_keep, ]
    num_rp <- nrow(kept_pairs)

    # Relabeling the records that are still candidates to match
    kept_records <- c(kept_pairs[, 1], kept_pairs[, 2])
    # not_included_records <- setdiff(1:sum(file_sizes), kept_records)
    labels <- cbind(kept_records, as.numeric(factor(kept_records)))

    kept_pairs[, 1] <- labels[1:num_rp, 2]
    kept_pairs[, 2] <- labels[(num_rp + 1):(2 * num_rp), 2]

    new_file_sizes <- rep(0, comparison_list$K)
    cum_file_sizes <- c(0, cumsum(file_sizes))

    for(k in 1:comparison_list$K){
        which_pairs <- kept_records > cum_file_sizes[k] &
            kept_records < (cum_file_sizes[k + 1] + 1)
        new_file_sizes[k] <- length(unique(labels[which_pairs, 1]))
    }


    # Calculate statistics needed for the Gibbs sampler
    new_file_labels <- c()
    for(k in 1:K){ new_file_labels <- c(new_file_labels,
                                        rep(k, new_file_sizes[k])) }
    num_fp <- K * (K + 1) / 2
    L <- sum(field_levels)

    # fp_matrix <- create_fp_matrix(K)
    temp_fp <- expand.grid(1:K, 1:K)[2:1]
    temp_fp <- temp_fp[temp_fp[, 1] <= temp_fp[, 2], ]
    fp_matrix <- matrix(0, nrow = K, ncol = K)
    for(i in 1:nrow(temp_fp)){
        fp_matrix[temp_fp[i, 1], temp_fp[i, 2]] <- i
        fp_matrix[temp_fp[i, 2], temp_fp[i, 1]] <- i
    }

    # rp_to_fp: a matrix that says which record pairs belong to which file pairs
    # rp_to_fp[fp, rp] is TRUE if the records being compared in comp_data[rp, ]
    # belong to the file pair fp, and is FALSE otherwise
    new_rp_to_fp <- matrix(FALSE, nrow = num_fp, ncol = num_rp)
    for(rp in 1:num_rp){
        i <- kept_pairs[rp, 1]
        j <- kept_pairs[rp, 2]
        new_rp_to_fp[fp_matrix[new_file_labels[i], new_file_labels[j]], rp] <- T
    }
    # ab: a vector of length L*num_fp which gives the count for each
    # disagreement level for each field for each file pair in comparisons
    new_ab <- rep(0, num_fp * L)
    ones <- rep(1, num_rp)
    for(fp in 1:num_fp){
        new_ab[((fp - 1) * L + 1):(fp * L)] <- ones[new_rp_to_fp[fp, ]] %*%
            kept_comparisons[new_rp_to_fp[fp, ], ]
    }

    list(record_pairs = kept_pairs, comparisons = kept_comparisons, K = K,
         file_sizes = new_file_sizes, duplicates = comparison_list$duplicates,
         field_levels = comparison_list$field_levels,
         file_labels = new_file_labels,  rp_to_fp = new_rp_to_fp, ab = new_ab,
         file_sizes_not_included = comparison_list$file_sizes - new_file_sizes,
         ab_not_included = comparison_list$ab - new_ab,
         labels = labels, pairs_to_keep = pairs_to_keep, cc = cc)

}
