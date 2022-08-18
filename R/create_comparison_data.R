#' Create Comparison Data
#'
#' Create comparison data for all pairs of records, except for those records in
#' files which are assumed to have no duplicates.
#'
#' The purpose of this function is to construct comparison vectors for each pair
#' of records. In order to construct these vectors, one needs to specify the
#' \code{types} and \code{breaks} arguments. The \code{types} argument specifies
#' how each field should be compared, and the \code{breaks} argument specifies
#' how to discretize these comparisons.
#'
#' Currently, the \code{types} argument supports three types of field
#' comparisons: binary, absolute difference, and the normalized Levenshtein
#' distance. Please contact the package maintainer if you need a new type of
#' comparison to be supported.
#'
#' The \code{breaks} argument should be a \code{list}, with with one element for
#' each field. If a field is being compared with a binary comparison, i.e.
#' \code{types[f]="bi"}, then the corresponding element of \code{breaks} should
#' be \code{NA}, i.e. \code{breaks[[f]]=NA}. If a field is being compared with a
#' numeric or string comparison, then the corresponding element of \code{breaks}
#' should be a vector of cut points used to discretize the comparisons. To give
#' more detail, suppose you pass in cut points
#' \code{breaks[[f]]=c(cut_1, ...,cut_L)}. These cut points
#' discretize the range of the comparisons into \code{L+1} intervals:
#' \eqn{I_0=(-\infty, cut_1], I_1=(cut_1, cut_2], ..., I_L=(cut_L, \infty]}. The
#' raw comparisons, which lie in \eqn{[0,\infty)} for numeric comparisons and
#' \eqn{[0,1]} for  string comparisons, are then replaced with indicators of
#' which interval the comparisons lie in. The interval \eqn{I_0} corresponds to
#' the lowest level of disagreement for a comparison, while the interval
#' \eqn{I_L} corresponds to the highest level of disagreement for a comparison.
#'
#' @param records A \code{data.frame} containing the records to be linked, where
#' each column of \code{records} is a field to be compared. If there are
#' multiple files, \code{records} should be obtained by stacking the files on
#' top of each other so that \code{records[1:file_sizes[1], ]} contains the
#' records for file \code{1},
#' \code{records[(file_sizes[1] + 1):(file_sizes[1] + file_sizes[2]), ]}
#' contains the records for file \code{2}, and so on. Missing values should be
#' coded as \code{NA}.
#' @param types A \code{character} vector, indicating the comparison to be used
#' for each field (i.e. each column of \code{records}). The options are:
#' \code{"bi"} for binary comparisons, \code{"nu"} for numeric comparisons
#' (absolute difference), \code{"lv"} for string comparisons (normalized
#' Levenshtein distance), \code{"lv_sep"} for string comparisons (normalized
#' Levenshtein distance) where each string may contain multiple spellings
#' separated by the "|" character. We assume that fields using options
#'  \code{"bi"}, \code{"lv"}, and \code{"lv_sep"} are  of class
#'  \code{character}, and fields using the \code{"nu"} option are of class
#'  \code{numeric}. For fields using the \code{"lv_sep"} option, for each record
#' pair the normalized Levenshtein distance is computed between each possible
#' spelling, and the minimum normalized Levenshtein distance between spellings
#' is then used as the comparison for that record pair.
#' @param breaks A \code{list}, the same length as \code{types}, indicating the
#' break points used to compute disagreement levels for each fields'
#' comparisons. If \code{types[f]="bi"}, \code{breaks[[f]]} is ignored (and thus
#' can be set to \code{NA}). See Details for more information on specifying this
#' argument.
#' @param file_sizes A \code{numeric} vector indicating the size of each file.
#' @param duplicates A \code{numeric} vector indicating which files are assumed
#' to have duplicates. \code{duplicates[k]} should be \code{1} if file \code{k}
#' has duplicates, and \code{duplicates[k]} should be \code{0} if file \code{k}
#' has no duplicates. If any files do not have duplicates, we strongly recommend
#' that the largest such file is organized to be the first file.
# @param K The number of files, assumed to be of class \code{numeric}. This
# should equal the length of file sizes. K = length(file_sizes)
#' @return a list containing:
#' \describe{
#'   \item{\code{record_pairs}}{A \code{data.frame}, where each row
#'   contains the pair of records being compared in the corresponding row of
#'   \code{comparisons}. The rows are sorted in ascending order according to the
#'   first column, with ties broken according to the second column in ascending
#'   order. For any given row, the first column is less than the second column,
#'   i.e. \code{record_pairs[i, 1] < record_pairs[i, 2]} for each row \code{i}.}
#'   \item{\code{comparisons}}{A \code{logical} matrix, where each row contains
#'   the comparisons for the record pair in the corresponding row of
#'   \code{record_pairs}. Comparisons are in the same order as the columns of
#'   \code{records}, and are represented by \code{L + 1} columns of
#'   \code{TRUE/FALSE} indicators, where \code{L + 1} is the number of
#'   disagreement levels for the field based on \code{breaks}.}
#'   \item{\code{K}}{The number of files, assumed to be of class
#'   \code{numeric}.}
#'   \item{\code{file_sizes}}{A \code{numeric} vector of length \code{K},
#'   indicating the size of each file.}
#'   \item{\code{duplicates}}{A \code{numeric} vector of length \code{K},
#'   indicating which files are assumed to have duplicates. \code{duplicates[k]}
#'   should be \code{1} if file \code{k} has duplicates, and
#'   \code{duplicates[k]} should be \code{0} if file \code{k} has no duplicates.
#'   If any files do not have duplicates, we strongly recommend that the largest
#'   such file is organized to be the first file.}
#'   \item{\code{field_levels}}{A \code{numeric} vector indicating the number of
#'   disagreement levels for each field.}
#'   \item{\code{file_labels}}{An \code{integer} vector of length
#'   \code{sum(file_sizes)}, where \code{file_labels[i]} indicates which file
#'   record \code{i} is in.}
#'   \item{\code{fp_matrix}}{An \code{integer} matrix, where
#'   \code{fp_matrix[k1, k2]} is a label for the file pair \code{(k1, k2)}. Note
#'   that \code{fp_matrix[k1, k2] = fp_matrix[k2, k1]}.}
#'   \item{\code{rp_to_fp}}{A \code{logical} matrix that indicates which record
#'   pairs belong to which file pairs. \code{rp_to_fp[fp, rp]} is \code{TRUE} if
#'   the records  \code{record_pairs[rp, ]} belong to the file pair \code{fp},
#'   and is FALSE otherwise. Note that \code{fp} is given by the labelling in
#'   \code{fp_matrix}.}
#'   \item{\code{ab}}{An \code{integer} vector, of length
#'   \code{ncol(comparisons) * K * (K + 1) / 2} that indicates how many record
#'   pairs there are with a given disagreement level for a given field, for each
#'   file pair.}
#'   \item{\code{file_sizes_not_included}}{A \code{numeric} vector of \code{0}s.
#'   This element is non-zero when \code{\link{reduce_comparison_data}} is
#'   used.}
#'   \item{\code{ab_not_included}}{A \code{numeric} vector of \code{0}s. This
#'   element is non-zero when \code{\link{reduce_comparison_data}} is used.}
#'   \item{\code{labels}}{\code{NA}. This element is not \code{NA} when
#'   \code{\link{reduce_comparison_data}} is used.}
#'   \item{\code{pairs_to_keep}}{\code{NA}. This element is not \code{NA} when
#'   \code{\link{reduce_comparison_data}} is used.}
#'   \item{\code{cc}}{\code{0}. This element is non-zero when
#'   \code{\link{reduce_comparison_data}} is used.}
#' }
#' @export
#' @examples
#' ## Example with no duplicate dataset
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
#' ## Example with duplicate dataset
#' data(dup_data)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(dup_data$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = dup_data$file_sizes,
#'  duplicates = c(1, 1, 1))
create_comparison_data <- function(records, types, breaks, file_sizes,
                                   duplicates){
    # Input Checks
    if(!is.data.frame(records)){ stop("'records' is not a data frame") }
    if(length(types) != ncol(records)){
        stop("'types' is not the same length as the number of fields
             'ncol(records)'")
    }
    if(length(types) != ncol(records)){
        stop("'types' is not the same length as the number of fields
             'ncol(records)'")
    }
    if(length(breaks) != ncol(records)){
        stop("'breaks' is not the same length as the number of fields
             'ncol(records)'")
    }
    if(nrow(records) != sum(file_sizes)){
        stop("Number of records 'nrow(records)' is not the same as
             'sum(file_sizes)'")
    }
    if(length(duplicates) != length(file_sizes)){
        stop("'length(duplicates)' is not the same as 'length(file_sizes)'")
    }
    if(!is.numeric(file_sizes)){
        stop("'file_sizes' must be numeric")
    }
    if(!is.numeric(duplicates)){
        stop("'duplicates' must be numeric")
    }
    for(k in 1:length(file_sizes)){
        if(duplicates[k] != 0 & duplicates[k] != 1){
            stop("Elements of 'duplicates' must be 0 or 1")
        }
    }
    for(f in 1:ncol(records)){
        if(types[f] == "bi"){
            if(!is.character(records[, f])){
                stop(paste0("Column ",  f, " of 'records' must be of type
                            character if a binary comparison is going to be
                            used"))
            }
        }
        else if (types[f] == "lv"){
            if(!is.character(records[, f])){
                stop(paste0("Column ",  f, " of 'records' must be of type
                            character if a string comparison is going to be
                            used"))
            }
            if(sum(is.na(breaks[[f]])) > 0){
                stop(paste0("Element ",  f, " of 'breaks' must be specified if a
                            string comparison is going to be used"))
            }
        }
        else if (types[f] == "lv_sep"){
            if(!is.character(records[, f])){
                stop(paste0("Column ",  f, " of 'records' must be of type
                            character if a string comparison is going to be
                            used"))
            }
            if(sum(is.na(breaks[[f]])) > 0){
                stop(paste0("Element ",  f, " of 'breaks' must be specified if a
                            string comparison is going to be used"))
            }
        }
        else if (types[f] == "nu"){
            if(!is.numeric(records[, f])){
                stop(paste0("Column ",  f, " of 'records' must be of type
                            numeric if a numeric comparison is going to be
                            used"))
            }
            if(sum(is.na(breaks[[f]])) > 0){
                stop(paste0("Element ",  f, " of 'breaks' must be specified if a
                            numeric comparison is going to be used"))
            }
        }
        else if (types[f] == "longlat"){
            if(!is.character(records[, f])){
                stop(paste0("Column ",  f, " of 'records' must be of type
                            character if a comparison between longitudes and
                            latitudes is going to be used"))
            }
            if(sum(is.na(breaks[[f]])) > 0){
                stop(paste0("Element ",  f, " of 'breaks' must be specified if a
                            comparison between longitudes and  latitudes is
                            going to be used"))
            }
        }
        else{
            stop("Elements of 'types' must be one of 'bi', 'lv',  'lv_sep', or
                 'nu'")
        }
    }

    r <- nrow(records)
    FF <- ncol(records)
    K <- length(file_sizes)

    # Determine the valid record pairs according to duplicates
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
    colnames(record_pairs) <- c("i", "j")
    record_pairs <- as.data.frame(record_pairs)

    rps1 <- record_pairs[, 1]
    rps2 <- record_pairs[, 2]
    num_rp <- nrow(record_pairs)

    # Create the disagreement levels based on the comparisons
    comparisons <- matrix(0, nrow = num_rp, ncol = FF)
    colnames(comparisons) <- colnames(records)
    field_levels <- rep(0, FF)
    field_names <- colnames(records)
    names(field_levels) <- field_names
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
            raw_comp <- 1 - RecordLinkage::levenshteinSim(records[rps1, f],
                                                          records[rps2, f])
            # as.numeric(utils::adist(records[, f], records[, f]))/
            #     pmax(nchar(records[, f])[rps1], nchar(records[, f])[rps2])
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, f] <- cut(raw_comp, breaks = temp_breaks,
                                    labels = 1:(length(temp_breaks) - 1))
            field_levels[f] <- length(temp_breaks) - 1
        }
        if(types[f] == "lv_sep"){
            # Find how many spellings there are of field f in each record
            num_spellings <- stringr::str_count(records[, f],
                                                pattern =
                                                    stringr::fixed("|")) + 1
            # Find the maximum number of spellings of field f in a record
            max_spellings <- max(num_spellings)
            # Split field f into max_spellings columns
            split_records <-
                as.data.frame(stringr::str_split_fixed(records[, f],
                                                       pattern =
                                                           stringr::fixed("|"),
                                                       n = max_spellings),
                              stringsAsFactors = FALSE)
            split_records <- apply(split_records, 2, trimws)

            # Calculate Levenshtein distance between each spelling
            raw_comps <- matrix(0, nrow = num_rp,
                                ncol = choose(max_spellings, 2) + max_spellings)
            counter <- 0
            for(i in 1:max_spellings){
                for(j in i:max_spellings){
                    counter <- counter + 1
                    raw_comps[, counter] <- 1 -
                        RecordLinkage::levenshteinSim(split_records[rps1, i],
                                                      split_records[rps2, j])
                }
            }

            # Take minimum of Levenshtein distances for each record pair
            raw_comps[is.nan(raw_comps)] <- 1
            raw_comp <- rep(0, num_rp)
            for(i in 1:num_rp){
                raw_comp[i] <- min(raw_comps[i, ])
            }
            temp_breaks <- unique(c(-Inf, breaks[[f]], Inf))
            comparisons[, f] <- cut(raw_comp, breaks = temp_breaks,
                                    labels = 1:(length(temp_breaks) - 1))
            field_levels[f] <- length(temp_breaks) - 1
        }
        if(types[f] == "longlat"){
            split_records <-
                as.data.frame(stringr::str_split_fixed(records[, f],
                                                       pattern = ",", n = 2),
                              stringsAsFactors = FALSE)
            split_records$V1 <- as.numeric(split_records$V1)
            split_records$V2 <- as.numeric(split_records$V2)
            # Convert distance from meters to miles
            raw_comp <- geosphere::distGeo(split_records[rps1, ],
                                           split_records[rps2, ]) * 0.000621371
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
    L <- sum(field_levels)
    cum_levels <- c(0, cumsum(field_levels))
    new_comparisons <- matrix(FALSE, nrow = num_rp, ncol = L)
    new_col_names <- rep("", L)
    for(f in 1: FF){
        for(l in 1:field_levels[f]){
            new_col_names[cum_levels[f] + l] <- paste0(field_names[f], "_DL_",
                                                       l - 1)
            new_comparisons[, cum_levels[f] + l] <- comparisons[, f] == l
        }
    }
    colnames(new_comparisons) <- new_col_names
    new_comparisons[is.na(new_comparisons)] <- FALSE
    comparisons <- new_comparisons


    # Calculate summaries needed for the Gibbs sampler
    file_labels <- c()
    for(k in 1:K){ file_labels <- c(file_labels, rep(k, file_sizes[k])) }
    num_fp <- K * (K + 1) / 2

    temp_fp <- expand.grid(1:K, 1:K)[2:1]
    temp_fp <- temp_fp[temp_fp[, 1] <= temp_fp[, 2], ]
    fp_matrix <- matrix(0, nrow = K, ncol = K)
    for(i in 1:nrow(temp_fp)){
        fp_matrix[temp_fp[i, 1], temp_fp[i, 2]] <- i
        fp_matrix[temp_fp[i, 2], temp_fp[i, 1]] <- i
    }

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
    ab_names <- rep("", num_fp * L)
    for(fp in 1:num_fp){
        ab_names[((fp - 1) * L + 1):(fp * L)] <-
            paste0(new_col_names, "_FP_(", temp_fp[fp, 1], ",", temp_fp[fp, 2],
                   ")")
        ab[((fp - 1) * L + 1):(fp * L)] <- ones[rp_to_fp[fp, ]] %*%
            comparisons[rp_to_fp[fp, ], ]
    }
    names(ab) <- ab_names
    ab_not_included <- rep(0, num_fp * L)
    names(ab_not_included) <- ab_names

    list(record_pairs = record_pairs, comparisons = comparisons, K = K,
         file_sizes = file_sizes, duplicates = duplicates,
         field_levels = field_levels, file_labels = file_labels,
         fp_matrix = fp_matrix, rp_to_fp = rp_to_fp, ab = ab,
         file_sizes_not_included = rep(0, K), ab_not_included = ab_not_included,
         labels = NA, pairs_to_keep = NA, cc = 0)
}
