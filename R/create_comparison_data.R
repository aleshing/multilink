#' Create Comparison Data
#'
#' Description goes here.
#'
#' @param records The records to be linked, of class \code{data.frame}.
#' It's assumed that each column of \code{records} is a field to be compared,
#' and that if there are multiple files (i.e. \code{K}>1), \code{records} is
#' obtained by stacking the files on top of each other so that
#' \code{records[1:file_sizes[1], ]} contains the records for file \code{1},
#' \code{records[(file_sizes[1] + 1):(file_sizes[1] + file_sizes[2]), ]}
#' contains the records for file \code{2}, and so on. Missing values should be
#' coded as \code{NA}.
#' @param types A \code{character} vector, indicating the comparison to be used
#' for each column of \code{records}. The options are: \code{"bi"} for
#' binary comparisons, \code{"nu"} for numeric comparisons
#' (absolute difference), and \code{"lv"} for string comparisons (normalized
#' Levenshtein distance). We assume that fields using options \code{"bi"} and
#' \code{"lv"} are  of class \code{character}, and fields using the \code{"nu"}
#' option are of class \code{numeric}.
#' @param breaks A \code{list}, the same length as \code{types}, indicating the
#' break points used to compute disagreement levels for each fields'
#' comparisons. If \code{types[f]=="bi"}, \code{breaks[[f]]} is ignored.
#' @param K The number of files, assumed to be of class \code{numeric}.
#' @param file_sizes A \code{numeric} vector of length \code{K}, indicating the
#' size of each file.
#' @param duplicates A \code{numeric} vector of length \code{K}, indicating
#' which files have duplicates. \code{duplicates[k]} should be \code{1} if file
#' \code{k} has duplicates, and \code{duplicates[k]} should be \code{0} if file
#' \code{k} has no duplicates. If any files do not have duplicates, we strongly
#' recommend that the largest such file is organized to be the first file.
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
#' @examples
create_comparison_data <- function(records, types, breaks, K, file_sizes,
                                   duplicates){
    r <- nrow(records)
    FF <- ncol(records)

    # Determine the valid record pairs based on duplicates
    if(sum(duplicates) == K){ record_pairs <- t(utils::combn(r, 2)) }
    else{
        record_pairs <- data.frame(Var2 = c(), Var1 = c())
        file_sizes_cum <- cumsum(c(0, file_sizes))
        for(k in 1:K){
            temp_records <- (file_sizes_cum[k] + 1):file_sizes_cum[k+1]
            if(duplicates[k] == 1){
                temp_record_pairs <- t(utils::combn(temp_records, 2))
                colnames(temp_record_pairs) <- c("Var2", "Var1")
                record_pairs <- rbind(record_pairs, temp_record_pairs)
            }
            if(k < K){
                temp_records_rest <- (file_sizes_cum[k+1] + 1):r
                record_pairs <- rbind(record_pairs,
                                      expand.grid(temp_records_rest,
                                                  temp_records)[, 2:1])
            }
        }
    }
    rps1 <- record_pairs[, 1]
    rps2 <- record_pairs[, 2]
    num_rp <- nrow(record_pairs)

    # Create the disagreement levels based on the comparisons
    comparisons <- matrix(0, nrow = num_rp, ncol = FF)
    colnames(comparisons) <- colnames(records)
    field_levels <- rep(0, FF)
    for(f in 1:FF){
        print(paste0("Creating comparison data for field ",
                     colnames(records)[f]))
        if(types[f] == "bi"){
            raw_comp <- records[rps1, f] != records[rps2, f]
            comparisons[, f] <- as.numeric(raw_comp) + 1
            field_levels[f] <- 2
        }
        if(types[f] == "nu"){
            raw_comp <- abs(records[rps1, f] - records[rps2, f])
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, f] <- cut(raw_comp, breaks = temp_breaks,
                                    labels = 1:(length(temp_breaks) - 1))
            field_levels[f] <- length(temp_breaks) - 1
        }
        if(types[f] == "lv"){
            raw_comp <- 1 - RecordLinkage::levenshteinSim(records[rps1, f], records[rps2, f])
            # as.numeric(utils::adist(records[, f], records[, f]))/
            #     pmax(nchar(records[, f])[rps1], nchar(records[, f])[rps2])
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, f] <- cut(raw_comp, breaks = temp_breaks,
                                    labels = 1:(length(temp_breaks) - 1))
            field_levels[f] <- length(temp_breaks) - 1
        }
    }
    comparisons <- as.data.frame(comparisons)

    # Turn the disagreement levels into an observation matrix, i.e. a logical
    # matrix with of size num_rp x L, where each row corresponds to the observed
    # values of each field for each pair of records being compared. In
    # particular obs_mat[rp, c] is TRUE if the records being compared in
    # comp_data[rp, ] observe the field at the level which c corresponds to.
    # Replacing NAs by FALSE is justified by MAR and CI assumptions.
    xfun <- function(f) paste("(comparisons[,", f, "]==", 1:field_levels[f],
                              ")")
    expr1 <- paste(sapply(sapply(1:FF, xfun), FUN = paste, collapse = ","),
                   collapse = ",")
    expr2 <- paste("cbind(", expr1, ")")
    comparisons <- eval(parse(text = expr2))
    comparisons[is.na(comparisons)] <- FALSE

    # Calculate statistics needed for the Gibbs sampler
    file_labels <- c()
    for(k in 1:K){ file_labels <- c(file_labels, rep(k, file_sizes[k])) }
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
    rp_to_fp <- matrix(FALSE, nrow = num_fp, ncol = num_rp)
    for(rp in 1:num_rp){
        i <- rps1[rp]
        j <- rps2[rp]
        rp_to_fp[fp_matrix[file_labels[i], file_labels[j]], rp] <- TRUE
    }
    # ab: a vector of length L*num_fp which gives the count for each
    # disagreement level for each field for each file pair in comparisons
    ab <- rep(0, num_fp * L)
    ones <- rep(1, num_rp)
    for(fp in 1:num_fp){
        ab[((fp - 1) * L + 1):(fp * L)] <- ones[rp_to_fp[fp, ]] %*%
            comparisons[rp_to_fp[fp, ], ]
    }

    list(record_pairs = record_pairs, comparisons = comparisons, K = K,
         file_sizes = file_sizes, duplicates = duplicates,
         field_levels = field_levels, file_labels = file_labels,
         rp_to_fp = rp_to_fp, ab = ab, file_sizes_not_included = rep(0, K),
         ab_not_included = rep(0, num_fp * L), labels = NA, pairs_to_keep = NA,
         cc = 0)
}
