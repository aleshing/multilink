#ifndef SAMPLE_Z_H
#define SAMPLE_Z_H

#include <RcppArmadillo.h>
using namespace Rcpp;

int bin_to_int_rcpp(const arma::uvec& bin, const arma::vec& powers);

List sample_Z_rcpp(arma::vec Z, arma::mat clust_sizes, int n, arma::vec cont,
                   const arma::vec& log_like, const arma::vec& alphas,
                   int alpha_0, const arma::vec& dup_upper_bound,
                   List dup_count_prior, const arma::vec& n_prior, int r,
                   int r_1, const arma::mat& valid_rp,
                   const arma::vec& singleton_ind, const arma::umat& rp_ind,
                   const arma::vec& file_labels, const arma::vec& powers,
                   int flat, int no_dups, int cc, arma::umat Z_members,
                   arma::vec clust_sizes_collapsed, int indexing_used);

#endif //SAMPLE_Z_H
