#ifndef SAMPLE_Z_H
#define SAMPLE_Z_H

#include <RcppArmadillo.h>
using namespace Rcpp;

int bin_to_int_rcpp(const arma::uvec& bin, const arma::vec& powers);

int uniform_discrete(int N);

IntegerVector uniform_chaperones(const arma::vec& file_size_cum,
                                 const arma::umat& valid_fp_matrix,
                                 const arma::vec& fp_probs);

IntegerVector nonuniform_chaperones(const arma::vec& file_size_cum,
                                    const arma::umat& valid_fp_matrix,
                                    const arma::vec& fp_probs,
                                    const arma::umat& comparisons_chap,
                                    int num_comps,
                                    int num_field,
                                    const arma::umat& record_pairs,
                                    const arma::vec& comparison_rps_probs);

List sample_Z_rcpp(arma::vec Z, arma::mat clust_sizes, int n, arma::vec cont,
                   const arma::vec& log_like, const arma::vec& alphas,
                   int alpha_0, const arma::vec& dup_upper_bound,
                   List dup_count_prior, const arma::vec& n_prior, int r,
                   int r_1, const arma::mat& valid_rp,
                   const arma::vec& singleton_ind, const arma::umat& rp_ind,
                   const arma::vec& file_labels, const arma::vec& powers,
                   int flat, int no_dups, int cc, arma::umat Z_members,
                   arma::vec clust_sizes_collapsed, int indexing_used);

List sample_Z_rcpp_chaperones(arma::vec Z, arma::mat clust_sizes, int n, arma::vec cont,
                              const arma::vec& log_like, const arma::vec& alphas,
                              int alpha_0, const arma::vec& dup_upper_bound,
                              List dup_count_prior, const arma::vec& n_prior, int r,
                              int r_1, const arma::mat& valid_rp,
                              const arma::vec& singleton_ind, const arma::umat& rp_ind,
                              const arma::vec& file_labels, const arma::vec& powers,
                              int flat, int no_dups, int cc, arma::umat Z_members,
                              arma::vec clust_sizes_collapsed, int indexing_used,
                              int num_chap_iter, int num_rp, int chap_type,
                              const arma::vec& file_size_cum, const arma::umat& valid_fp_matrix,
                              const arma::vec& fp_probs,
                              const arma::umat& record_pairs,
                              List comparison_rps,
                              int comparison_rps_length, int extra_gibbs,
                              int num_restrict,
                              const arma::umat& comparisons_chap,
                              int num_field, const arma::vec& comparison_rps_probs);

#endif //SAMPLE_Z_H
