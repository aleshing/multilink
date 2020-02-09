#include "compute_log_like.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This funciton can be shoved into sample_phi
arma::vec compute_log_like_rcpp(arma::vec m, arma::vec u, int num_rp,
                                int num_fp, int L, arma::mat rp_to_fp,
                                arma::mat obs_mat) {
    // Is there a better way to construct this matrix?
    // And to subsequently compute the dot products we want?
    arma::vec log_ratio = log(m / u);
    // log_ratio_mat: Same size as obs_mat, where each row has the
    // appropriate log ratio for the pair of records in the corresponding
    // row in obs_mat
    arma::mat log_ratio_mat(num_rp, L);
    for(int fp = 0; fp < num_fp; fp++){
        arma::uvec ids = arma::find(rp_to_fp.row(fp));
        arma::rowvec temp = log_ratio.subvec(fp * L, (fp +1 ) * L - 1).t();
        for(int i = 0; i < ids.size(); i++){
            log_ratio_mat.row(ids(i)) = temp;
        }
    }
    arma::vec log_like = (obs_mat % log_ratio_mat) * arma::ones(L);
    return(log_like);
}
