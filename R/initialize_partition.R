#' Initialize the Partition
#'
#' Description goes here.
#'
#' @param comparison_list the output from a call to
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}
#' @param pairs_to_keep
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
initialize_partition <- function(comparison_list, pairs_to_keep, seed=NA){

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
                                    prob=rep(1, files_counts[j]))
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
