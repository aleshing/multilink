#' Relabel the Bayes Estimate of a Partition
#'
#' Relabel the Bayes estimate of a partition, for use after using indexing to
#' reduce the number of record pairs that are potential matches.
#'
#' When the function \code{\link{reduce_comparison_data}} is used to reduce the
#' number of record pairs that are potential matches, it may be the case that
#' some records are declared to not be potential matches to any other records.
#' In this case, the indexing method has made the decision that these records
#' have no matches, and thus we can remove them from the data set and relabel
#' the remaining records; see the documentation for \code{labels} in
#' \code{\link{reduce_comparison_data}} for information on how to go between the
#' original labeling and the new labeling. The purpose of this function is to
#' relabel the output of \code{\link{find_bayes_estimate}} when the function
#' \code{\link{reduce_comparison_data}} is used, so that the user doesn't have
#' to do this relabeling themselves.
#'
#' @param reduced_comparison_list The output from a call to
#' \code{\link{reduce_comparison_data}}.
#' @param bayes_estimate The output from a call to
#' \code{\link{find_bayes_estimate}}.
#'
#' @return
#' A \code{data.frame}, with as many rows as
#' \code{sum(reduced_comparison_list$file_sizes +
#' reduced_comparison_list$file_sizes_not_included)}, i.e. the number of
#' records originally input to \code{\link{create_comparison_data}}, before
#' indexing occurred. This \code{data.frame} has two columns,
#' \code{"original_labels"} and \code{"link_id"}. Given row \code{i} of
#' \code{records} originally input to \code{\link{create_comparison_data}},
#' the linkage id according to \code{bayes_estimate} is given by the \code{i}th
#' row of the \code{link_id} column. See the documentation for
#' \code{\link{find_bayes_estimate}} for information on how to interpret this
#' linkage id.
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\doi{https://doi.org/10.1080/01621459.2021.2013242}][\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#' @export
#'
#' @examples
#' # Example with small duplicate dataset
#' data(dup_data_small)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(dup_data_small$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = dup_data_small$file_sizes,
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
#' \dontrun{
#' results <- gibbs_sampler(reduced_comparison_list, prior_list, n_iter = 1000,
#'  seed = 42)
#'
#' # Find the full Bayes estimate
#' full_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = Inf, max_cc_size = 50)
#'
#' # Find the partial Bayes estimate
#' partial_estimate <- find_bayes_estimate(results$partitions, burn_in = 100,
#'  L_FNM = 1, L_FM1 = 1, L_FM2 = 2, L_A = 0.1, max_cc_size = 12)
#'
#' # Relabel the full and partial Bayes estimates
#' full_estimate_relabel <- relabel_bayes_estimate(reduced_comparison_list,
#'  full_estimate)
#'
#' partial_estimate_relabel <- relabel_bayes_estimate(reduced_comparison_list,
#'  partial_estimate)
#'
#' # Add columns to the records corresponding to their full and partial
#' # Bayes estimates
#' dup_data$records <- cbind(dup_data$records,
#'  full_estimate_id = full_estimate_relabel$link_id,
#'  partial_estimate_id = partial_estimate_relabel$link_id)
#' }
relabel_bayes_estimate <- function(reduced_comparison_list, bayes_estimate){
    # Input checks
    num_records_kept <- sum(reduced_comparison_list$file_sizes)
    num_records_not_included <-
        sum(reduced_comparison_list$file_sizes_not_included)
    length_estimate <- length(bayes_estimate)
    if(num_records_kept != num_records_kept){
        stop("Length of 'bayes_estimate' does not equal the number of records
             kept in 'reduced_comparison_list'. Make sure
             'reduced_comparison_list' is the same as was input to the function
             gibbs_sampler.")
    }

    label_dict <- as.data.frame(reduced_comparison_list$labels)
    link_estimates <- data.frame(new_labels = 1:length(bayes_estimate),
                                 link_id = bayes_estimate)
    num_records <- num_records_kept + num_records_not_included
    temp <- merge(reduced_comparison_list$labels, link_estimates,
                  by = "new_labels")
    temp$original_labels <- temp$kept_records
    temp <- temp[, c("original_labels", "link_id")]
    if(sum(reduced_comparison_list$file_sizes_not_included) != 0){
        temp <- rbind(temp,
                      data.frame(link_id = NA,
                                 original_labels =
                                     setdiff(1:num_records,
                                             label_dict$kept_records)))


        temp$link_id[(num_records_kept + 1):num_records] <-
            (max(bayes_estimate) + 1):(max(bayes_estimate) +
                                           (num_records - num_records_kept))
    }
    temp <- temp[order(temp$original_labels),]

    return(temp)
}
