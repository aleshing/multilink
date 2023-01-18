#' Initialize the Partition
#'
#' Generate an initialization for the partition in the case when it is assumed
#' there are no duplicates in all files (so that the partition is a matching).
#'
#' When it is assumed that there are no duplicates in all files, and
#' \code{\link{reduce_comparison_data}} is not used to reduce the number of
#' potential matches, the Gibbs sampler used for posterior inference may
#' experience slow mixing when using an initialization for the partition where
#' each record is in its own cluster (the default option for the Gibbs sampler).
#' The purpose of this function is to provide an alternative initialization
#' scheme.
#'
#' To use this initialization scheme, the user passes in a \code{logical} vector
#' that indicates which record pairs are potential matches according to an
#' indexing method (as in \code{\link{reduce_comparison_data}}). Note that this
#' indexing is only used to generate the initialization, it is not used for
#' inference. The initialization scheme first finds the transitive closure of
#' the potential matches, which partitions the records into blocks. Within each
#' block of records, the scheme randomly selects a record from each file, and
#' these selected records are then placed in the same cluster for the partition
#' initialization. All other records are placed in their own clusters.
#'
#' @param comparison_list the output from a call to
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}.
#' Note that in order to correctly specify the initialization, if
#' \code{\link{reduce_comparison_data}} is used to the reduce the number of
#' record pairs that are candidate matches, then the output of
#' \code{\link{reduce_comparison_data}} (not
#' \code{\link{create_comparison_data}}) should be used for this argument.
#' @param pairs_to_keep A \code{logical} vector, the same length as
#' \code{comparison_list$record_pairs}, indicating which record pairs are
#' potential matches in the initialization.
#' @param seed The seed to use to generate the initialization.
#'
#' @return an \code{integer} vector of arbitrary labels of length
#' \code{sum(comparison_list$file_sizes)}, giving an initialization for the
#' partition.
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\href{https://doi.org/10.1080/01621459.2021.2013242}{Published}] [\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#'
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
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = no_dup_data$file_sizes,
#'  duplicates = c(0, 0, 0))
#'
#' # Find initialization for the matching
#' # Only keep pairs of records as potential matches in the initialization for
#' # which neither gname nor fname disagree at the highest level
#' pairs_to_keep <- (comparison_list$comparisons[, "gname_DL_3"] != TRUE) &
#'  (comparison_list$comparisons[, "fname_DL_3"] != TRUE)
#' Z_init <- initialize_partition(comparison_list, pairs_to_keep, seed = 42)
initialize_partition <- function(comparison_list, pairs_to_keep, seed = NA){
    # Input Checks
    if(!is.logical(pairs_to_keep)){
        stop("'pairs_to_keep' must be logical")
    }
    if(length(pairs_to_keep) != nrow(comparison_list$record_pairs)){
        stop("'pairs_to_keep' is not the same length as the number of record
             pairs 'nrow(comparison_list$record_pairs)'")
    }

    record_pairs <- comparison_list$record_pairs
    file_sizes <- comparison_list$file_sizes

    pairs_to_keep[is.na(pairs_to_keep)] <- FALSE
    pairs <- record_pairs[pairs_to_keep, ]
    g <- igraph::graph_from_edgelist(as.matrix(pairs), directed=FALSE)
    components <- igraph::components(g)

    Z <- components$membership
    free <- max(Z) + 1
    # If the tail of the records is made up of singletons, need to add them to Z
    if(length(Z) != sum(file_sizes)){
        Z[(max(pairs) + 1):sum(file_sizes)] <-
            free:(free - 1 + sum(file_sizes) - max(pairs))

    }
    free <- max(Z) + 1

    K <- comparison_list$K
    file_labels <- comparison_list$file_labels

    for(i in 1:components$no){
        if(components$csize[i] > 1){
            members <- which(components$membership == i)
            files <- file_labels[members]
            files_counts <- rep(0, K)
            for(j in 1:K){
                files_counts[j] <- sum(files == j)
            }
            old_members <- c()
            for(j in 1:K){
                if(files_counts[j] > 1){
                    k <- 1
                    if(!is.na(seed)){
                        set.seed(seed)
                        k <- sample(1:files_counts[j], size = 1, replace = TRUE,
                                    prob = rep(1, files_counts[j]))
                    }
                    old_members <- c(old_members, members[files == j][-k])
                }
            }
            Z[old_members] <- -1
        }
    }
    Z[Z == -1] <- free:(free + sum(Z == - 1) - 1)
    return(Z)
}
