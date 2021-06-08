#include "sample_phi.h"
#include "sample_Z.h"
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List gibbs_loop_rcpp(int n_iter, arma::mat Z_samp, arma::mat clust_sizes_samp,
                     arma::mat cont_samp, arma::mat m_samp, arma::mat u_samp,
                     const arma::vec& mus, const arma::vec& nus,
                     const arma::vec& alphas, int alpha_0,
                     const arma::vec& dup_upper_bound, List dup_count_prior,
                     const arma::vec& n_prior, arma::vec cont,
                     arma::mat clust_sizes, int n, const arma::vec& ab,
                     const arma::mat& obs_mat, const arma::umat& record_pairs,
                     int flat, int r, int r_1, const arma::mat& valid_rp,
                     const arma::vec& singleton_ind, const arma::umat& rp_ind,
                     const arma::vec& file_labels, const arma::vec& powers,
                     int L, int num_fp, int num_rp, int num_field,
                     const arma::mat& rp_to_fp, const arma::vec& level_cum,
                     int no_dups, const arma::vec& valid_fp, int cc,
                     arma::umat Z_members, arma::vec clust_sizes_collapsed,
                     int indexing_used, int single_likelihood,
                     const arma::vec& single_nus, const arma::vec& single_ab){
    for(int i = 1; i < n_iter; i++){
        //// Every 10 iterations check to see if the user has interrupted the
        //// sampler
        if(i % 10 == 0 || i == 1){
            checkUserInterrupt();
        }
        if(i % 100 == 0 || i == 1){
            Rcout << "Beginning iteration " << i << "/" << n_iter - 1 <<
                std::endl;
        }

        //// Inner Gibbs loop: Fixed scan of comparison data model parameters
        // coref_vec: a binary vector where a 1 represents a record pair being
        // in the same cluster in the current partition
        arma::vec Z_curr = Z_samp.col(i - 1);
        arma::uvec coref_vec_u = Z_curr.elem(record_pairs.col(0) - 1) ==
            Z_curr.elem(record_pairs.col(1) - 1);
        arma::vec coref_vec = arma::conv_to<arma::vec>::from(coref_vec_u);

        List phi = sample_phi_rcpp(coref_vec, obs_mat, ab, mus, nus, L,
                                   num_fp, num_rp, num_field, rp_to_fp,
                                   level_cum, valid_fp, single_likelihood,
                                   single_nus, single_ab);

        m_samp.col(i) = as<arma::vec>(phi["m"]);
        u_samp.col(i) = as<arma::vec>(phi["u"]);

        // log_like: the log likelihood contribution for each record pair
        arma::vec log_like = as<arma::vec>(phi["log_like"]);

        //// Inner Gibbs loop: Fixed scan of cluster assignments for records
        //// Use labels from the end of iteration i-1 for the start of iteration
        //// i
        List Z = sample_Z_rcpp(Z_curr, clust_sizes, n, cont, log_like, alphas,
                               alpha_0, dup_upper_bound, dup_count_prior,
                               n_prior, r, r_1, valid_rp, singleton_ind, rp_ind,
                               file_labels, powers, flat, no_dups, cc,
                               Z_members, clust_sizes_collapsed, indexing_used);

        Z_samp.col(i) = as<arma::vec>(Z["Z"]);
        cont = as<arma::vec>(Z["cont"]);
        n = as<int>(Z["n"]);
        clust_sizes = as<arma::mat>(Z["clust_sizes"]);
        cont_samp.col(i) = cont;
        Z_members = as<arma::umat>(Z["Z_members"]);
        clust_sizes_collapsed = as<arma::vec>(Z["clust_sizes_collapsed"]);
        clust_sizes_samp.col(i) = clust_sizes_collapsed;
    }

    return(List::create(Named("m_samp") = m_samp, Named("u_samp") = u_samp,
                        Named("Z_samp") = Z_samp,
                        Named("cont_samp") = cont_samp,
                        Named("clust_sizes_samp") = clust_sizes_samp));
}
