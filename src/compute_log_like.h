#ifndef COMPUTE_LOG_LIKE_H
#define COMPUTE_LOG_LIKE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec compute_log_like_rcpp(arma::vec m, arma::vec u, int num_rp,
                                int num_fp, int L, arma::mat rp_to_fp,
                                arma::mat obs_mat);

#endif //COMPUTE_LOG_LIKE_H
