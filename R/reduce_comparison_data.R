#' Reduce Comparison Data Size
#'
#' Use indexing to reduce the number of record pairs that are potential matches.
#'
#' When using comparison-based record linkage methods, scalability is a concern,
#' as the number of record pairs is quadratic in the number of records. In
#' order to address these concerns, it's common to declare certain record pairs
#' to not be potential matches a priori, using indexing methods. The user is
#' free to index using any method they like, as long as they can produce a
#' \code{logical} vector that indicates which record pairs are potential matches
#' according to their indexing method. We recommend, if the user chosen indexing
#' method does not output potential matches that are transitive, to set the
#' \code{cc} argument to \code{1}. By transitive we mean, for any three records
#' \eqn{i}, \eqn{j}, and \eqn{k}, if \eqn{i} and \eqn{j} are potential matches,
#' and \eqn{j} and \eqn{k} are potential matches, then \eqn{i} and \eqn{k} are
#' potential matches. Non-transitive indexing schemes can lead to poor mixing of
#' the Gibbs sampler used for posterior inference, and suggests that the
#' indexing method used may have been too stringent.
#'
#' If indexing is used, it may be the case that some records are declared to not
#' be potential matches to any other records. In this case, the indexing method
#' has made the decision that these records have no matches, and thus we can
#' remove them from the data set and relabel the remaining records; see the
#' documentation for \code{labels} for information on how to go between the
#' original labeling and the new labeling.
#'
#' If indexing is used, comparisons for record pairs that aren't potential
#' matches are still used during inference, where they're used to inform the
#' distribution of comparisons for non-matches.
#'
#' @param comparison_list The output of a call to
#' \code{\link{create_comparison_data}}.
#' @param pairs_to_keep A \code{logical} vector, the same length as
#' \code{comparison_list$record_pairs}, indicating which record pairs should be
#' kept as potential matches. These potential matches do not have to be
#' transitive (see the argument \code{cc}).
#' @param cc A \code{numeric} indicator of whether to find the transitive
#' closure of \code{pairs_to_keep}, and use these potential matches instead
#' of just those from \code{pairs_to_keep}. \code{cc} should be \code{1} if the
#' transitive closure is being used, and \code{cc} should be \code{0} if the
#' transitive closure is not being used. We recommend setting \code{cc} to
#' \code{1}.
#' @return a list containing:
#' \describe{
#'   \item{\code{record_pairs}}{A \code{data.frame}, where each row
#'   contains the pair of records being compared in the corresponding row of
#'   \code{comparisons}. The rows are sorted in ascending order according to the
#'   first column, with ties broken according to the second column in ascending
#'   order. For any given row, the first column is less than the second column,
#'   i.e. \code{record_pairs[i, 1] < record_pairs[i, 2]} for each row \code{i}.
#'   If according to \code{pairs_to_keep} there are records which are not
#'   potential matches to any other records, the remaining records are
#'   relabeled (see \code{labels}).}
#'   \item{\code{comparisons}}{A \code{logical} matrix, where each row contains
#'   the comparisons between the record pair in the corresponding row of
#'   \code{record_pairs}. Comparisons are in the same order as the columns of
#'   \code{records}, and are represented by \code{L + 1} columns of
#'   \code{TRUE/FALSE} indicators, where \code{L + 1} is the number of
#'   disagreement levels for the field based on \code{breaks}.}
#'   \item{\code{K}}{The number of files, assumed to be of class
#'   \code{numeric}.}
#'   \item{\code{file_sizes}}{A \code{numeric} vector of length \code{K},
#'   indicating the size of each file. If according to \code{pairs_to_keep}
#'   there are records which are not potential matches to any other records, the
#'   remaining records are relabeled (see \code{labels}), and \code{file_sizes}
#'   now represents the sizes of each file after removing such records.}
#'   \item{\code{duplicates}}{A \code{numeric} vector of length \code{K},
#'   indicating which files are assumed to have duplicates. \code{duplicates[k]}
#'   should be \code{1} if file \code{k} has duplicates, and
#'   \code{duplicates[k]} should be \code{0} if file \code{k} has no
#'   duplicates.}
#'   \item{\code{field_levels}}{A \code{numeric} vector indicating the number of
#'   disagreement levels for each field.}
#'   \item{\code{file_labels}}{An \code{integer} vector of length
#'   \code{sum(file_sizes)}, where \code{file_labels[i]} indicated which file
#'   record \code{i} is in.}
#'   \item{\code{fp_matrix}}{An \code{integer} matrix, where
#'   \code{fp_matrix[k1, k2]} is a label for the file pair \code{(k1, k2)}. Note
#'   that \code{fp_matrix[k1, k2] = fp_matrix[k2, k1]}.}
#'   \item{\code{rp_to_fp}}{A \code{logical} matrix that indicates which record
#'   pairs belong to which file pairs. \code{rp_to_fp[fp, rp]} is \code{TRUE} if
#'   the records  \code{record_pairs[rp, ]} belong to the file pair \code{fp},
#'   and is FALSE otherwise. Note that \code{fp} is given by the labeling in
#'   \code{fp_matrix}.}
#'   \item{\code{ab}}{An \code{integer} vector, of length
#'   \code{ncol(comparisons) * K * (K + 1) / 2} that indicates how many record
#'   pairs there are with a given disagreement level for a given field, for each
#'   file pair.}
#'   \item{\code{file_sizes_not_included}}{If according to \code{pairs_to_keep}
#'   there are records which are not potential matches to any other records, the
#'   remaining records are relabeled (see \code{labels}), and
#'   \code{file_sizes_not_included} indicates, for each file, the number of such
#'   records that were removed.}
#'   \item{\code{ab_not_included}}{For record pairs not included according to
#'   \code{pairs_to_keep}, this is an \code{integer} vector, of length
#'   \code{ncol(comparisons) * K * (K + 1) / 2} that indicates how many record
#'   pairs there are with a given disagreement level for a given field, for each
#'   file pair.}
#'   \item{\code{labels}}{If according to \code{pairs_to_keep}
#'   there are records which are not potential matches to any other records, the
#'   remaining records are relabeled. \code{labels} provides a dictionary that
#'   indicates, for each of the new labels, which record in the original
#'   labeling the new label corresponds to. In particular, the first column
#'   indicates the record in the original labeling, and the second column
#'   indicates the new labeling.}
#'   \item{\code{pairs_to_keep}}{A \code{logical} vector, the same length as
#'   \code{comparison_list$record_pairs}, indicating which record pairs were
#'   kept as potential matches. This may not be the same as the input
#'   \code{pairs_to_keep} if \code{cc} was set to 1.}
#'   \item{\code{cc}}{A \code{numeric} indicator of whether the connected
#'   components of the potential matches are closed under transitivity.}
#' }
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\href{https://doi.org/10.1080/01621459.2021.2013242}{Published}] [\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#' @export
#'
#' @examples
#' # Example with duplicate dataset
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
#' # Use indexing to only keep pairs of records for which neither gname nor
#' # fname disagree at the highest level
#' pairs_to_keep <- (comparison_list$comparisons[, "gname_DL_3"] != TRUE) &
#'  (comparison_list$comparisons[, "fname_DL_3"] != TRUE)
#' reduced_comparison_list <- reduce_comparison_data(comparison_list,
#'  pairs_to_keep, cc = 1)
reduce_comparison_data <- function(comparison_list, pairs_to_keep, cc = 1){
    # Input Checks
    if(!is.logical(pairs_to_keep)){
        stop("'pairs_to_keep' must be logical")
    }
    if(length(pairs_to_keep) != nrow(comparison_list$record_pairs)){
        stop("'pairs_to_keep' is not the same length as the number of record
             pairs 'nrow(comparison_list$record_pairs)'")
    }
    if(!is.numeric(cc)){
        stop("'cc' must be numeric")
    }
    if(cc != 0 & cc!= 1){
        stop("'cc' must be 0 or 1")
    }


    # Dump the contents of comparison_list
    record_pairs <- comparison_list$record_pairs
    comparisons <- comparison_list$comparisons
    file_sizes <- comparison_list$file_sizes
    K <- comparison_list$K
    field_levels <- comparison_list$field_levels

    # Find the transitive closure of pairs_to_keep
    if(cc == 1){
        g <- igraph::graph_from_edgelist(
            as.matrix(record_pairs[pairs_to_keep, ]), directed = FALSE)
        components <- igraph::components(g)

        new_kept_pairs <- data.frame("i" = NULL, "j" = NULL)
        for(i in 1:length(components$csize)){
            if(components$csize[i] > 1){
                ind <- which(components$membership == i)
                # Note we don't need to construct these connected components
                # carefully (the next code chunk deals with this)
                new_kept_pairs <- rbind(new_kept_pairs, t(utils::combn(ind, 2)))
            }
        }

        colnames(new_kept_pairs) <- c("i", "j")
        colnames(record_pairs) <- c("i", "j")

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
    # This next line is actual magic, it preserves the ordering in the labels
    labels <- cbind(kept_records, new_labels = as.numeric(factor(kept_records)))

    kept_pairs[, 1] <- labels[1:num_rp, 2]
    kept_pairs[, 2] <- labels[(num_rp + 1):(2 * num_rp), 2]

    new_file_sizes <- rep(0, comparison_list$K)
    cum_file_sizes <- c(0, cumsum(file_sizes))

    for(k in 1:comparison_list$K){
        which_pairs <- kept_records > cum_file_sizes[k] &
            kept_records < (cum_file_sizes[k + 1] + 1)
        new_file_sizes[k] <- length(unique(labels[which_pairs, 1]))
    }

    # Calculate summaries needed for the Gibbs sampler
    new_file_labels <- c()
    for(k in 1:K){ new_file_labels <- c(new_file_labels,
                                        rep(k, new_file_sizes[k])) }
    num_fp <- K * (K + 1) / 2
    L <- sum(field_levels)

    temp_fp <- expand.grid(1:K, 1:K)[2:1]
    temp_fp <- temp_fp[temp_fp[, 1] <= temp_fp[, 2], ]
    fp_matrix <- matrix(0, nrow = K, ncol = K)
    for(i in 1:nrow(temp_fp)){
        fp_matrix[temp_fp[i, 1], temp_fp[i, 2]] <- i
        fp_matrix[temp_fp[i, 2], temp_fp[i, 1]] <- i
    }

    new_rp_to_fp <- matrix(FALSE, nrow = num_fp, ncol = num_rp)
    for(rp in 1:num_rp){
        i <- kept_pairs[rp, 1]
        j <- kept_pairs[rp, 2]
        new_rp_to_fp[fp_matrix[new_file_labels[i], new_file_labels[j]], rp] <- T
    }

    new_ab <- rep(0, num_fp * L)
    ones <- rep(1, num_rp)
    for(fp in 1:num_fp){
        new_ab[((fp - 1) * L + 1):(fp * L)] <- ones[new_rp_to_fp[fp, ]] %*%
            kept_comparisons[new_rp_to_fp[fp, ], ]
    }
    names(new_ab) <- names(comparison_list$ab)

    list(record_pairs = kept_pairs, comparisons = kept_comparisons, K = K,
         file_sizes = new_file_sizes, duplicates = comparison_list$duplicates,
         field_levels = comparison_list$field_levels,
         file_labels = new_file_labels,  fp_matrix = fp_matrix,
         rp_to_fp = new_rp_to_fp, ab = new_ab,
         file_sizes_not_included = comparison_list$file_sizes - new_file_sizes,
         ab_not_included = comparison_list$ab - new_ab,
         labels = labels[!duplicated(labels[, 1]), ],
         pairs_to_keep = pairs_to_keep, cc = cc)

}
