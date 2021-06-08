#include "sample_phi.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Draw a sample from a Dirichlet w/ parameter vector alpha
arma::vec rdir_rcpp(const arma::vec& alpha) {
    int size = alpha.n_elem;
    arma::vec sample = arma::zeros(size);
    double normalizer = 0;

    for(int i = 0; i < size; i++) {
        double temp = R::rgamma(alpha(i), 1.0);
        sample(i) = temp;
        normalizer += temp;
    }
    sample = sample / normalizer;

    return(sample);
}

List sample_phi_rcpp(const arma::vec& coref_vec, const arma::mat& obs_mat,
                     const arma::vec& ab, const arma::vec& mus,
                     const arma::vec& nus, int L, int num_fp, int num_rp,
                     int num_field, const arma::mat& rp_to_fp,
                     const arma::vec& level_cum, const arma::vec& valid_fp,
                     int single_likelihood, const arma::vec& single_nus,
                     const arma::vec& single_ab){
    arma::vec a = arma::zeros<arma::vec>(L);
    arma::vec m = arma::zeros<arma::vec>(num_fp * L);
    arma::vec u = arma::zeros<arma::vec>(num_fp * L);

    // Sample m and u while computing the log likelihood
    // log_like: The log likelihood for all record pairs
    arma::vec log_like(num_rp);

    if(single_likelihood == 0){
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
                    rdir_rcpp(a.subvec(start - fp * L, end - fp * L) +
                    mus.subvec(start, end));
                u.subvec(start, end) =
                    rdir_rcpp(ab.subvec(start, end) - a.subvec(start - fp * L,
                                        end - fp * L) + nus.subvec(start, end));
            }

            // Compute the log likelihood for records belonging to the current
            // file pair
            arma::rowvec log_ratio =
                log(m.subvec(fp * L, (fp + 1) * L - 1) /
                    u.subvec(fp * L, (fp + 1) * L - 1)).t();
            for(int i = 0; i < ids.size(); i++){
                log_like(ids(i)) = arma::sum(log_ratio % obs_mat.row(ids(i)));
            }
        }
    }
    else{
        // Compute a based on current partition
        a = obs_mat.t() * coref_vec;

        // Sample m and u
        for(int f = 0; f < num_field; f++){
            int start = level_cum(f) - 1;
            int end = level_cum(f + 1) - 2;
            m.subvec(start, end) = rdir_rcpp(a.subvec(start, end)
                                                 + mus.subvec(start, end));
            u.subvec(start, end) = rdir_rcpp(single_ab.subvec(start, end) -
            a.subvec(start, end) + single_nus.subvec(start, end));
        }

        // Compute the log likelihood
        arma::vec log_ratio =
            log(m.subvec(0, L - 1) / u.subvec(0, L - 1));
        log_like = obs_mat * log_ratio;
    }



    return List::create(Named("m") = m,
                        Named("u") = u,
                        Named("log_like") = log_like);
}
