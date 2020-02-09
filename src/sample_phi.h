#ifndef SAMPLE_PHI_H
#define SAMPLE_PHI_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat rdir_rcpp(int n, const arma::vec& alpha);

List sample_phi_rcpp(const arma::vec& coref_vec, const arma::mat& obs_mat,
                     const arma::vec& ab, const arma::vec& mus,
                     const arma::vec& nus, int L, int num_fp, int num_rp,
                     int num_field, const arma::mat& rp_to_fp,
                     const arma::vec& level_cum, const arma::vec& valid_fp);

#endif //SAMPLE_PHI_H
