#include "sample_phi.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Draw n samples from a Dirichlet w/ parameter vector alpha
// Cribbed from http://www.mjdenny.com/blog.html
arma::mat rdir_rcpp(int n, const arma::vec& alpha) {
    int size = alpha.n_elem;
    // Each row will be a draw from a Dirichlet
    arma::mat samps = arma::zeros(n, size);
    // Each pass of this loop is a draw from a Dirichlet
    for (int i=0; i<n; ++i) {
        double sum_term = 0;
        // Loop through the distribution and draw Gamma variables
        for (int j=0; j<size; ++j) {
            double cur = R::rgamma(alpha[j], 1.0);
            samps(i, j) = cur;
            sum_term += cur;
        }
        // Normalize
        for (int j=0; j<size; ++j) {
            samps(i, j) = samps(i, j)/sum_term;
        }
    }
    // Actually let's make each column a draw
    return(samps.t());
}

List sample_phi_rcpp(const arma::vec& coref_vec, const arma::mat& obs_mat,
                     const arma::vec& ab, const arma::vec& mus,
                     const arma::vec& nus, int L, int num_fp, int num_rp,
                     int num_field, const arma::mat& rp_to_fp,
                     const arma::vec& level_cum, const arma::vec& valid_fp){
    arma::vec a = arma::zeros<arma::vec>(L);
    arma::vec m = arma::zeros<arma::vec>(num_fp * L);
    arma::vec u = arma::zeros<arma::vec>(num_fp * L);

    // Sample m and u while computing the log likelihood
    // log_ratio_mat: The log likelihood for all record pairs
    arma::vec log_like(num_rp);
    for(int fp_ind = 0; fp_ind < valid_fp.size(); fp_ind++){
        int fp = valid_fp(fp_ind);

        // Compute a based on current partition for the current file pair
        arma::uvec ids = arma::find(rp_to_fp.row(fp));
        a = obs_mat.rows(ids).t() * coref_vec.elem(ids);

        // Sample m and u for the current file pair
        for(int f = 0; f < num_field; f++){
            int start = fp * L + level_cum(f) - 1;
            int end = fp * L + level_cum(f + 1) - 2;
            m.subvec(start, end) =
                rdir_rcpp(1, a.subvec(start - fp * L, end - fp * L) +
                mus.subvec(start, end));
            u.subvec(start, end) =
                rdir_rcpp(1, ab.subvec(start, end) - a.subvec(start - fp * L,
                                       end - fp * L) + nus.subvec(start, end));
        }

        // Compute the log likelihood for records belonging to the current file pair
        arma::rowvec log_ratio =
            log(m.subvec(fp * L, (fp + 1) * L - 1) /
            u.subvec(fp * L, (fp + 1) * L - 1)).t();
        for(int i = 0; i < ids.size(); i++){
            log_like(ids(i)) = arma::sum(log_ratio % obs_mat.row(ids(i)));
        }
    }


    return List::create(Named("m") = m,
                        Named("u") = u,
                        Named("log_like") = log_like);
}
