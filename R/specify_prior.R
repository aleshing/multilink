#' Specify the Prior Distributions
#'
#' Specify the prior distributions for the \eqn{m} and \eqn{u} parameters of the
#' models for comparison data among matches and non-matches, and the partition.
#'
#' The purpose of this function is to specify prior distributions for all
#' parameters of the model. Please note that if
#' \code{\link{reduce_comparison_data}} is used to the reduce the number of
#' record pairs that are potential matches, then the output of
#' \code{\link{reduce_comparison_data}} (not
#' \code{\link{create_comparison_data}}) should be used as input.
#'
#' For the hyperparameters of the Dirichlet priors for the \eqn{m}
#' and \eqn{u} parameters for the comparisons among matches and non-matches,
#' respectively, we recommend using a flat prior. This is accomplished by
#' setting \code{mus=NA} and \code{nus=NA}. Informative prior specifications
#' are possible, but in practice they will be overwhelmed by the large number of
#' comparisons.
#'
#' For the prior for partitions, we do not recommend using a flat prior. Instead
#' we recommend using our structure prior for partitions. By setting
#' \code{flat=0} and the remaining arguments to \code{NA}, one obtains the
#' default specification for the structured prior that we have found to perform
#' well in simulation studies. The structured prior for partitions is specified
#' as follows:
#' \itemize{
#'   \item Specify a prior for \code{n}, the number of clusters represented in
#'   the records. Note that this includes records determined not to be potential
#'   matches to any other records using \code{\link{reduce_comparison_data}}.
#'   Currently, a uniform prior and a scale prior for \code{n} are supported.
#'   Our default specification uses a uniform prior.
#'   \item Specify a prior for the overlap table (see the documentation for
#'   \code{alphas} for more information).  Currently a Dirichlet-multinomial
#'   prior is supported. Our default specification sets all hyperparameters of
#'   the Dirichlet-multinomial prior to \code{1}.
#'   \item For each file, specify a prior for the number of duplicates in each
#'   cluster. As a part of this prior, we specify the maximum number of records
#'   in a cluster for each file, through \code{dup_upper_bound}. When there
#'   are assumed to be no duplicates in a file, the maximum number of records in
#'   a cluster for that file is set to \code{1}. When there are assumed to be
#'   duplicates in a file, we recommend setting the maximum number of records in
#'   a cluster for that file to be less than the file size, if prior knowledge
#'   allows. Currently, a Poisson prior for the the number of duplicates in
#'   each cluster is supported. Our default specification uses a Poisson prior
#'   with mean \code{1}.
#' }
#' Please contact the package maintainer if you need new prior families
#' for \code{n} or the number of duplicates in each cluster to be supported.
#'
#' @param comparison_list the output from a call to
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}.
#' Note that in order to correctly specify the prior, if
#' \code{\link{reduce_comparison_data}} is used to the reduce the number of
#' record pairs that are potential matches, then the output of
#' \code{\link{reduce_comparison_data}} (not
#' \code{\link{create_comparison_data}}) should be used for this argument.
#' @param mus,nus The hyperparameters of the Dirichlet priors for the \eqn{m}
#' and \eqn{u} parameters for the comparisons among matches and non-matches,
#' respectively. These are positive \code{numeric} vectors which have length
#' equal to the number of  columns of \code{comparison_list$comparisons} times
#' the number of file pairs
#' \code{(comparison_list$K * (comparison_list$K + 1) / 2)}. If set to
#' \code{NA}, flat priors are used. We recommend using flat priors for \eqn{m}
#' and \eqn{u}.
#' @param flat A \code{numeric} indicator of whether a flat prior for partitions
#' should be used.  \code{flat} should be \code{1} if a flat prior is used, and
#' \code{flat} should be \code{0} if a structured prior is used. If a flat prior
#' is used, the remaining arguments should be set to \code{NA}. Otherwise, the
#' remaining arguments should be specified. We do not recommend using a flat
#' prior for partitions in general.
#' @param alphas The hyperparameters for the Dirichlet-multinomial overlap table
#' prior, a positive \code{numeric} vector of length
#' \code{2 ^ comparison_list$K - 1}. The indexing of these hyperparameters is
#' based on the the \code{comparison_list$K}-bit binary representation of the
#' inclusion patterns of the overlap table. To give a few examples, suppose
#' \code{comparison_list$K} is \code{3}. \code{1} in \code{3}-bit binary is
#' \code{001}, so \code{alphas[1]} is the hyperparameter for the
#' \code{001} cell of the overlap table, representing clusters containing only
#' records from the third file. \code{2} in \code{3}-bit binary is
#' \code{010}, so \code{alphas[2]} is the hyperparameter for the
#' \code{010} cell of the overlap table, representing clusters containing only
#' records from the second file. \code{3} in \code{3}-bit binary is
#' \code{011}, so \code{alphas[3]} is the hyperparameter for the
#' \code{011} cell of the overlap table, representing clusters containing only
#' records from the second and third files. If set to \code{NA}, the
#' hyperparameters will all be set to \code{1}.
#' @param dup_upper_bound A \code{numeric} vector indicating the maximum number
#' of duplicates, from each file, allowed in each cluster. For a given file
#' \code{k}, \code{dup_upper_bound[k]} should be between \code{1} and
#' \code{comparison_list$file_sizes[k]}, i.e. even if you don't want to impose
#' an upper bound, you have to implicitly place an upper bound: the number of
#' records in a file. If set to \code{NA}, the upper bound for file \code{k}
#' will be set to \code{1} if no duplicates are allowed for that file, or
#' \code{comparison_list$file_sizes[k]} if duplicates are allowed for that file.
#' @param dup_count_prior_family A \code{character} vector indicating the
#' prior distribution family used for the number of duplicates in each cluster,
#' for each file. Currently the only option is \code{"Poisson"} for a Poisson
#' prior, truncated to lie between \code{1} and \code{dup_upper_bound[k]}. The
#' mean parameter of the Poisson distribution is specified using the
#' \code{dup_count_prior_pars} argument. If set to \code{NA}, a Poisson prior
#' with mean \code{1} will be used.
#' @param dup_count_prior_pars A \code{list} containing the parameters for
#' the prior distribution for the number of duplicates in each cluster, for each
#' file. For file \code{k}, when \code{dup_count_prior_family[k]="Poisson"},
#' \code{dup_count_prior_pars[[k]]} is a positive constant representing the mean
#' of the Poisson prior.
#' @param n_prior_family A \code{character} indicating the prior distribution
#' family used for \code{n}, the number of clusters represented in the
#' records. Note that this includes records determined not to be potential
#' matches to any other records using \code{\link{reduce_comparison_data}}.
#' Currently the there are two options: \code{"uniform"} for a uniform prior
#' for \code{n}, i.e. \eqn{p(n) \propto 1}, and \code{"scale"} for a scale prior
#' for \code{n}, i.e. \eqn{p(n) \propto 1/n}. If set to \code{NA}, a uniform
#' prior will be used.
#' @param n_prior_pars Currently set to \code{NA}. When more prior distribution
#' families for \code{n} are implemented, this will be a vector of parameters
#' for those priors.
#' @return a list containing:
#' \describe{
#'   \item{\code{mus}}{The hyperparameters of the Dirichlet priors for the
#'   \code{m} parameters for the comparisons among matches.}
#'   \item{\code{nus}}{The hyperparameters of the Dirichlet priors for the
#'   \code{u} parameters for the comparisons among non-matches. Includes data
#'   from comparisons of record pairs that were declared to not be potential
#'   matches using \code{\link{reduce_comparison_data}}.}
#'   \item{\code{flat}}{A \code{numeric} indicator of whether a flat prior for
#'   partitions should be used. \code{flat} is \code{1} if a flat prior is used,
#'   and \code{flat} is \code{0} if a structured prior is used.}
#'   \item{\code{no_dups}}{A \code{numeric} indicator of whether no duplicates
#'   are allowed in all of the files.}
#'   \item{\code{alphas}}{The hyperparameters for the Dirichlet-multinomial
#'   overlap table prior, a positive \code{numeric} vector of length
#'   \code{2 ^ comparison_list$K}, where the first element is \code{0}.}
#'   \item{\code{alpha_0}}{The sum of \code{alphas}.}
#'   \item{\code{dup_upper_bound}}{A \code{numeric} vector indicating the
#'   maximum number of duplicates, from each file, allowed in each cluster. For
#'   a given file \code{k}, \code{dup_upper_bound[k]} should be between \code{1}
#'   and \code{comparison_list$file_sizes[k]}, i.e. even if you don't want to
#'   impose an upper bound, you have to implicitly place an upper bound: the
#'   number of records in a file.}
#'   \item{\code{log_dup_count_prior}}{A \code{list} containing the log density
#'   of the prior distribution for the number of duplicates in each cluster, for
#'   each file.}
#'   \item{\code{log_n_prior}}{A \code{numeric} vector containing the log
#'   density of the prior distribution for the number of clusters represented in
#'   the records.}
#'   \item{\code{nus_specified}}{The \code{nus} before data from comparisons of
#'   record pairs that were declared to not be potential matches using
#'   \code{\link{reduce_comparison_data}} are added. Used for input checking.}
#' }
#' @references Serge Aleshin-Guendel & Mauricio Sadinle (2022). Multifile Partitioning for Record Linkage and Duplicate Detection. \emph{Journal of the
#' American Statistical Association}. [\doi{https://doi.org/10.1080/01621459.2021.2013242}] [\href{https://arxiv.org/abs/2110.03839}{arXiv}]
#' @export
#'
#' @examples
#' # Example with small no duplicate dataset
#' data(no_dup_data_small)
#'
#' # Create the comparison data
#' comparison_list <- create_comparison_data(no_dup_data_small$records,
#'  types = c("bi", "lv", "lv", "lv", "lv", "bi", "bi"),
#'  breaks = list(NA,  c(0, 0.25, 0.5),  c(0, 0.25, 0.5),
#'                c(0, 0.25, 0.5), c(0, 0.25, 0.5),  NA, NA),
#'  file_sizes = no_dup_data_small$file_sizes,
#'  duplicates = c(0, 0, 0))
#'
#' # Specify the prior
#' prior_list <- specify_prior(comparison_list, mus = NA, nus = NA, flat = 0,
#'  alphas = rep(1, 7), dup_upper_bound = c(1, 1, 1),
#'  dup_count_prior_family = NA, dup_count_prior_pars = NA,
#'  n_prior_family = "uniform", n_prior_pars = NA)
#'
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
specify_prior <- function(comparison_list, mus = NA, nus = NA, flat = 0,
                          alphas = NA, dup_upper_bound = NA,
                          dup_count_prior_family = NA,
                          dup_count_prior_pars = NA,
                          n_prior_family = NA, n_prior_pars = NA){
    # Input checks
    if(!is.na(mus)){
        if(!is.numeric(mus)){
            stop("'mus' must be numeric")
        }
        if(length(mus) != ncol(comparison_list$comparisons)){
            stop("''length(mus) and 'ncol(comparison_list$comparisons)' must be
                 the same")
        }
    }
    if(!is.na(nus)){
        if(!is.numeric(nus)){
            stop("'nus' must be numeric")
        }
        if(length(nus) != ncol(comparison_list$comparisons)){
            stop("''length(nus) and 'ncol(comparison_list$comparisons)' must be
                 the same")
        }
    }
    if(!is.numeric(flat)){
        stop("'flat' must be numeric")
    }
    if(flat != 0 & flat!= 1){
        stop("'flat' must be 0 or 1")
    }
    if(sum(!is.na(alphas)) > 0){
        if(!is.numeric(alphas)){
            stop("'alphas' must be numeric")
        }
        if(length(alphas) != 2 ^ comparison_list$K - 1){
            stop("''length(alphas) and '2 ^ comparison_list$K - 1' must be the
                 same")
        }
    }
    if(sum(!is.na(dup_upper_bound)) > 0){
        if(!is.numeric(dup_upper_bound)){
            stop("'dup_upper_bound' must be numeric")
        }
        for(k in 1:comparison_list$K){
            if(dup_upper_bound[k] < 1 | dup_upper_bound[k] >
               comparison_list$file_sizes[k]){
                stop(paste0("Element", k,
                " of 'dup_upper_bound' must be between 1 and element ", k,
                " of 'comparison_list$file_sizes'"))
            }
            if(dup_upper_bound[k] > 1 & comparison_list$duplicates[k] == 0){
                stop(paste0("Element", k,
                " of 'dup_upper_bound' must be 1 if element ", k,
                " of 'comparison_list$duplicates' is 0"))
            }
        }
    }
    for(k in 1:comparison_list$K){
        if(!is.na(dup_count_prior_family[k])){
            if(dup_count_prior_family[k] != "Poisson"){
                stop(paste0("Element", k, " of 'dup_count_prior_family' must be
                            'Poisson' or NA currently"))
            }
            if(is.na(dup_count_prior_pars[[k]])){
                stop(paste0("Element", k, " of 'dup_count_prior_pars' must be
                            specified"))
            }
            if(!is.numeric(dup_count_prior_pars[[k]])){
                stop(paste0("Element", k, " of 'dup_count_prior_pars' must be
                            numeric"))
            }
            if(length(dup_count_prior_pars[[k]]) != 1){
                stop(paste0("Element", k, " of 'dup_count_prior_pars' must be
                            length 1"))
            }
        }
    }
    if(!is.na(n_prior_family)){
        if(n_prior_family != "uniform" & n_prior_family != "scale"){
            stop("'n_prior_family' must be 'uniform', 'scale', or NA currently")
        }
    }


    K <- comparison_list$K
    num_fp <- K * (K + 1) / 2
    L <- sum(comparison_list$field_levels)

    # Specify flat priors for m and/or u if mus and/or nus were set to NA
    if(sum(is.na(mus)) > 0){ mus <- rep(1, num_fp * L)}
    if(sum(is.na(nus)) > 0){ nus <- rep(1, num_fp * L)}

    nus_specified <- nus
    #### This next line means you need to respecify the prior anytime you change
    #### the comparison data
    #### Add in prior non-coreferent counts using ab_full
    nus <- nus + comparison_list$ab_not_included
    names(mus) <- names(comparison_list$ab)
    names(nus_specified) <- names(comparison_list$ab)

    if(flat){
        upper <- comparison_list$file_sizes
        for(k in 1:K){
            if(comparison_list$duplicates[k] == 0){
                upper[k] <- 1
            }
        }
        return(list(mus = mus, nus = nus, flat = flat,
                    no_dups = as.numeric(sum(comparison_list$duplicates) == 0),
                    alphas = NA, alpha_0 = NA,
                    dup_upper_bound = upper,
                    log_dup_count_prior = NA, log_n_prior = NA,
                    nus_specified = nus_specified))
    }
    else{
        # Specify the prior for n
        r <- sum(comparison_list$file_sizes) +
            sum(comparison_list$file_sizes_not_included)
        log_n_prior <- rep(0, r)
        if(is.na(n_prior_family)){ n_prior_family <- "uniform" }
        if(n_prior_family == "uniform"){
            log_n_prior <- rep(-log(r), r)
        }
        else if(n_prior_family == "scale"){
            log_n_prior <- -log(1:r) - log(sum(1 / (1:r)))
        }

        # Specify the prior for the overlap table
        if(sum(is.na(alphas)) > 0){
            alphas <- rep(1, (2 ^ comparison_list$K) - 1)
        }
        alphas <- c(0, alphas)
        X <- matrix(as.numeric(intToBits(0:(2 ^ comparison_list$K - 1))),
                    nrow = 2 ^ comparison_list$K,
                    byrow = TRUE)[, comparison_list$K:1]
        if(K == 1){
            names(alphas) <- c("0", "1")
        }
        else{
            names(alphas) <- apply(X, 1, paste0, collapse = "")
        }


        # If the maximum number of duplicates for each file were set to NA,
        # specify implicit upper bounds
        if(sum(!is.na(dup_upper_bound)) == 0){
            upper <- comparison_list$file_sizes
            for(k in 1:K){
                if(comparison_list$duplicates[k] == 0){
                    upper[k] <- 1
                }
            }
            dup_upper_bound <- upper
        }

        # Specify prior for duplicate counts
        log_dup_count_prior <- vector("list", comparison_list$K)
        for(k in 1:K){
            if(!comparison_list$duplicates[k]){ log_dup_count_prior[[k]] <- 0 }
            else{
                if(is.na(dup_count_prior_family[k])){
                    dup_count_prior_family[k] <- "Poisson"
                    dup_count_prior_pars[[k]] <- 1
                }
                if(dup_count_prior_family[k] == "Poisson"){
                    temp_prior <- stats::dpois(1:dup_upper_bound[k],
                                               dup_count_prior_pars[[k]])
                    log_dup_count_prior[[k]] <- log(temp_prior /
                                                        sum(temp_prior))
                }

            }
        }

        return(list(mus = mus, nus = nus, flat = flat,
                    no_dups = as.numeric(sum(comparison_list$duplicates) == 0),
                    alphas = alphas, alpha_0 = sum(alphas),
                    dup_upper_bound = dup_upper_bound,
                    log_dup_count_prior = log_dup_count_prior,
                    log_n_prior = log_n_prior, nus_specified = nus_specified))
    }
}
