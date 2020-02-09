#' Specify the Partition Prior
#'
#' Description goes here.
#'
#' @param comparison_list the output from a call to
#' \code{\link{create_comparison_data}} or \code{\link{reduce_comparison_data}}
#' @param flat
#' @param alphas
#' @param dup_upper_bound
#' @param dup_count_prior_family
#' @param dup_count_prior_pars
#' @param n_prior_family
#' @param n_prior_pars
#'
#' @return a list containing:
#' \describe{
#'   \item{\code{flat}}{}
#'   \item{\code{no_dups}}{}
#'   \item{\code{alphas}}{}
#'   \item{\code{alpha_0}}{}
#'   \item{\code{dup_upper_bound}}{}
#'   \item{\code{log_dup_count_prior}}{}
#'   \item{\code{log_n_prior}}{}
#' }
#' @export
#'
#' @examples
specify_partition_prior <- function(comparison_list, flat,
                                    alphas = NA, dup_upper_bound = NA,
                                    dup_count_prior_family = NA,
                                    dup_count_prior_pars = NA,
                                    n_prior_family = NA, n_prior_pars = NA){
    if(flat){
        upper <- comparison_list$file_sizes
        for(k in 1:comparison_list$K){
            if(comparison_list$duplicates[k] == 0){
                upper[k] <- 1
            }
        }
        return(list(flat = flat,
                    no_dups = as.numeric(sum(comparison_list$duplicates) == 0),
                    alphas = NA, alpha_0 = NA,
                    dup_upper_bound = upper,
                    log_dup_count_prior = NA, log_n_prior = NA))
    }
    else{
        r <- sum(comparison_list$file_sizes) +
            sum(comparison_list$file_sizes_not_included)
        log_n_prior <- rep(0, r)
        if(is.na(n_prior_family)){ n_prior_family <- "uniform" }
        if(n_prior_family == "uniform"){
            log_n_prior <- log(rep(1, r) / r)
        }

        if(sum(is.na(alphas))>0){
            alphas <- rep(1, (2 ^ comparison_list$K) - 1)
        }
        alphas <- c(0, alphas)

        log_dup_count_prior <- vector("list", comparison_list$K)
        for(k in 1:comparison_list$K){
            if(!comparison_list$duplicates[k]){ log_dup_count_prior[[k]] <- 0 }
            else{
                if(is.na(dup_count_prior_family[k])){
                    dup_count_prior_family[k] <- "Poisson"
                    dup_count_prior_pars[[k]] <- c(1)
                }
                if(dup_count_prior_family[k] == "Poisson"){
                    temp_prior <- stats::dpois(1:dup_upper_bound[k],
                                               dup_count_prior_pars[[k]])
                    log_dup_count_prior[[k]] <- log(temp_prior /
                                                        sum(temp_prior))
                }

            }
        }

        return(list(flat = flat,
                    no_dups = as.numeric(sum(comparison_list$duplicates) == 0),
                    alphas = alphas, alpha_0 = sum(alphas),
                    dup_upper_bound = dup_upper_bound,
                    log_dup_count_prior = log_dup_count_prior,
                    log_n_prior = log_n_prior))
    }
}
