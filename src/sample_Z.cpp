#include "sample_Z.h"
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

int bin_to_int_rcpp(const arma::uvec& bin, const arma::vec& powers){
    arma::vec temp = bin.t()*powers;
    return(temp(0));
}

int uniform_discrete(int N){
    return(ceil(runif(1)(0) * N));
}

IntegerVector uniform_chaperones(const arma::vec& file_size_cum,
                                 const arma::umat& valid_fp_matrix,
                                 const arma::vec& fp_probs){
    IntegerVector chaperones = {0, 0};

    IntegerVector fps = seq(0, fp_probs.size() - 1);
    int fp_samp = RcppArmadillo::sample(fps, 1, TRUE, fp_probs)(0);

    int file_1 = valid_fp_matrix(0, fp_samp);
    int file_2 = valid_fp_matrix(1, fp_samp);
    int file_size_1 = file_size_cum(file_1) - file_size_cum(file_1 - 1);

    if(file_1 == file_2){
        chaperones(0) = uniform_discrete(file_size_1) + file_size_cum(file_1 - 1) - 1;
        chaperones(1) = uniform_discrete(file_size_1 - 1) + file_size_cum(file_1 - 1) - 1;
        if(chaperones(1) >= chaperones(0)){
            chaperones(1)++;
        }

    }
    else{
        int file_size_2 = file_size_cum(file_2) - file_size_cum(file_2 - 1);
        chaperones(0) = uniform_discrete(file_size_1) + file_size_cum(file_1 - 1) - 1;
        chaperones(1) = uniform_discrete(file_size_2 - 1) + file_size_cum(file_2 - 1) - 1;
    }

    return(chaperones);
}

IntegerVector nonuniform_chaperones(const arma::vec& file_size_cum,
                                    const arma::umat& valid_fp_matrix,
                                    const arma::vec& fp_probs,
                                    const arma::umat& record_pairs,
                                    List comparison_rps,
                                    int comparison_rps_length,
                                    const arma::vec& comparison_rps_probs){
    IntegerVector chaperones = {0, 0};
    IntegerVector rps = seq(-1, comparison_rps_length - 1);
    int opt = RcppArmadillo::sample(rps, 1, TRUE, comparison_rps_probs)(0);

    if(opt == -1){
        chaperones = uniform_chaperones(file_size_cum, valid_fp_matrix, fp_probs);
    }
    else{
        arma::vec candidate_rps = as<arma::vec>(comparison_rps[opt]);
        int rp = uniform_discrete(candidate_rps.size()) - 1;
        chaperones = record_pairs.row(candidate_rps[rp] - 1);
    }

    return(chaperones);
}

List sample_Z_rcpp(arma::vec Z, arma::mat clust_sizes, int n, arma::vec cont,
                   const arma::vec& log_like, const arma::vec& alphas,
                   int alpha_0, const arma::vec& dup_upper_bound,
                   List dup_count_prior, const arma::vec& n_prior, int r,
                   int r_1, const arma::mat& valid_rp,
                   const arma::vec& singleton_ind, const arma::umat& rp_ind,
                   const arma::vec& file_labels, const arma::vec& powers,
                   int flat, int no_dups, int cc, arma::umat Z_members,
                   arma::vec clust_sizes_collapsed, int indexing_used){
    for(int j = r_1; j < r; j++){
        // k: the file number of the current record
        int k = file_labels(j);

        // Get the label of the cluster record j is currently in
        int label = Z(j);

        // Update cont, clust_sizes, and n to reflect partition without
        // current record
        //
        // If record j is the only record in its cluster from file k, update
        // inclusion patterns
        if(clust_sizes(k - 1, label - 1) == 1){
            // Find the inclusion pattern of record j's cluster
            arma::uvec ind = clust_sizes.col(label - 1) > 0;
            int oldind = bin_to_int_rcpp(ind, powers);

            // Reduce the number of clusters with the inclusion pattern of
            // the cluster including record j
            cont(oldind)--;

            // If j wasn't in a singleton cluster, increase the number of
            // clusters with the inclusion pattern of the cluster without
            // record j
            if(sum(ind) > 1){
                ind(k - 1) = 0;
                int newind = bin_to_int_rcpp(ind, powers);
                cont(newind)++;
            }
            // Update n if j was in a singleton cluster

            else{ n--; }
        }

        // Using Z_members slows down inference when we have duplicates
        if(no_dups){
            int label_size = clust_sizes_collapsed(label - 1);
            // If j was in a singleton cluster, update Z_members
            if(label_size == 1){
                Z_members.submat(label - 1, 0, label - 1, 0) = 0;
            }
            // If j wasn't in a singleton cluster, update Z_members
            else if(label_size > 1){
                arma::urowvec temp_1 = Z_members.submat(label - 1, 0,
                                                        label - 1,
                                                        label_size - 1);
                temp_1.subvec(0, label_size - 2) =
                    temp_1.elem(arma::find(temp_1 != (j + 1))).t();
                temp_1.subvec(label_size - 1, label_size - 1) = 0;
                Z_members.submat(label - 1, 0,
                                 label - 1, label_size - 1) = temp_1;
            }
        }

        // Update clust_sizes
        clust_sizes(k - 1, label - 1)--;
        clust_sizes_collapsed(label - 1)--;

        // Give current record an empty label
        Z(j) = -1;
        arma::vec valid_c;
        arma::urowvec not_empty_clust;
        if(no_dups){
            // Valid clusters are just clusters where there are already records,
            // where none of the records are from file k, along with the empty
            // cluster
            not_empty_clust = clust_sizes_collapsed.t() > 0;
            if(indexing_used){
                if(cc){
                    arma::vec possible_c =
                        arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)))));

                    int empty_counter = 0;
                    for(int i = 0; i < possible_c.size(); i++){
                        if(possible_c(i) != -1){
                            if(clust_sizes(k - 1, possible_c(i) - 1) > 0){
                                possible_c(i) = -1;
                                empty_counter++;
                            }
                        }
                    }
                    valid_c = arma::unique(possible_c);
                }
                else{
                    arma::urowvec ind = clust_sizes.row(k - 1) == 0 && not_empty_clust;
                    valid_c = arma::conv_to<arma::vec>::from(arma::find(ind.t()) + 1);
                    arma::vec not_possible_c =
                        arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)
                                                                      == 0))));

                    int empty_counter = 0;
                    int npc_counter = 0;
                    int npc_size = not_possible_c.size();
                    for(int i = 0; i < valid_c.size(); i++){
                        int temp = valid_c(i);
                        if(npc_counter < npc_size){
                            int temp_npc = not_possible_c(npc_counter);
                            while(npc_counter < npc_size & temp_npc <= temp){
                                if(temp == temp_npc){
                                    valid_c(i) = -1;
                                    empty_counter++;
                                }
                                npc_counter++;
                                if(npc_counter < npc_size){
                                    temp_npc = not_possible_c(npc_counter);
                                }
                            }
                        }
                    }
                    if(empty_counter == 0){
                        int len = valid_c.size();
                        valid_c.resize(len + 1);
                        valid_c(len) = -1;
                    }
                    else{ valid_c = arma::unique(valid_c); }
                }
            }
            else{
                arma::urowvec ind = clust_sizes.row(k - 1) == 0 && not_empty_clust;
                valid_c = arma::conv_to<arma::vec>::from(arma::find(ind.t()) + 1);
                int len = valid_c.size();
                valid_c.resize(len + 1);
                valid_c(len) = -1;
            }
        }
        else{
            // Find valid clusters, where -1 from the empty label will
            // correspond to forming a new cluster. Make sure to sort them!
            arma::vec possible_c =
                arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)))));
            arma::vec not_possible_c;
            // Have to do this next line since some of the clusters w/ valid
            // records may also contain non-valid records if we don't have
            // transitive closures of all connected components. Make sure to
            // sort them!
            if(!cc){
                not_possible_c =
                    arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)
                                                                  == 0))));
            }

            // Can't join clusters that have more than the allowed number of
            // records from the same file! And make sure all of the clusters
            // in possible_c are actually possible if we don't have transitive
            // closures of all connected components
            if(!cc){
                int npc_counter = 0;
                int npc_size = not_possible_c.size();
                for(int i = 0; i < possible_c.size(); i++){
                    int temp = possible_c(i);
                    if(temp != -1){
                        if(clust_sizes(k - 1, temp - 1) >=
                           dup_upper_bound(k - 1)){
                            possible_c(i) = -1;
                        }
                        else if(npc_counter < npc_size){
                            int temp_npc = not_possible_c(npc_counter);
                            while(npc_counter < npc_size & temp_npc <= temp){
                                if(temp == temp_npc){
                                    possible_c(i) = -1;
                                }
                                npc_counter++;
                                if(npc_counter < npc_size){
                                    temp_npc = not_possible_c(npc_counter);
                                }
                            }
                        }
                    }
                }
            }
            else{
                for(int i = 0; i < possible_c.size(); i++){
                    int temp = possible_c(i);
                    if(temp != -1){
                        if(clust_sizes(k - 1, temp - 1) >=
                           dup_upper_bound(k - 1)){
                            possible_c(i) = -1;
                        }
                    }
                }
            }
            valid_c = arma::unique(possible_c);
        }

        // Calculate cluster assignment probabilities
        arma::vec log_probs = arma::zeros(valid_c.size());
        for(int c_index = 0; c_index < valid_c.size(); c_index++){
            // c_label: the label of the current cluster we're calculating
            // the log_probs for
            int c_label = valid_c(c_index);
            // If c_label isn't an empty cluster
            if(c_label != -1){
                // Get the records in c_label
                arma::uvec records;
                // Using Z_members slows down inference when we have duplicates
                if(no_dups){
                    records = (Z_members.submat(c_label - 1, 0, c_label - 1,
                                                clust_sizes_collapsed(c_label
                                                                          - 1)
                                                    - 1) - 1).t();
                }
                else{
                    records = arma::find(Z == c_label);
                }

                // Add the likelihood contribution
                arma::uvec pairs = rp_ind.col(j);
                pairs = pairs.elem(records) - 1;
                if(sum(pairs < 0) > 0){
                    Rcout << "Something went wrong, not getting the correct\
                    record pairs" << j << std::endl;
                }
                log_probs(c_index) += arma::sum(log_like.elem(pairs));

                if(!flat){
                    // c_k: the number of records from file k in the current
                    // cluster
                    int c_k = clust_sizes(k - 1, c_label - 1);
                    if(c_k >= dup_upper_bound(k - 1)){
                        Rcout << "Something went wrong, trying to sample a\
                        cluster with too many records from file: " << k <<
                            std::endl;
                    }
                    // If there are no other records from file k in c_label
                    if(c_k == 0){
                        arma::uvec ind = clust_sizes.col(c_label - 1) > 0;
                        int ind_without = bin_to_int_rcpp(ind, powers);
                        ind(k - 1) = 1;
                        int ind_with = bin_to_int_rcpp(ind, powers);

                        log_probs(c_index) +=
                            as<arma::vec>(dup_count_prior[k - 1])(0) +
                            std::log(cont(ind_with) + alphas(ind_with)) -
                            std::log(cont(ind_without) + alphas(ind_without)
                                         - 1);
                    }
                    // If there are other records from file k in c_label
                    else{
                        log_probs(c_index) +=
                            std::log(c_k + 1) +
                            as<arma::vec>(dup_count_prior[k - 1])(c_k) -
                            as<arma::vec>(dup_count_prior[k - 1])(c_k - 1);
                    }
                }
            }
            // If c_label is an empty cluster
            else{
                if(!flat){
                    log_probs(c_index) +=
                        as<arma::vec>(dup_count_prior(k - 1))(0) +
                        std::log(n + 1) +
                        std::log(cont(singleton_ind(k - 1) - 1) +
                        alphas(singleton_ind(k - 1) - 1)) -
                        std::log(n + alpha_0) + n_prior(n) - n_prior(n - 1);
                }
            }
        }

        //Sample the cluster assignments
        arma::vec probs = exp(log_probs);

        int samp = RcppArmadillo::sample(valid_c, 1, TRUE, probs)(0);
        if(samp == -1){
            if(no_dups){
                arma::uvec empty_clusts =
                    arma::find((not_empty_clust == 0)) + 1;
                samp = empty_clusts(0);
            }
            else{
                for(int find = 1; find < r+1; find++){
                    arma::uvec ind = Z == find;
                    if(sum(ind) == 0){
                        samp = find;
                        find = r+1;
                    }
                }
            }
        }
        Z(j) = samp;
        // Update cont, clust_sizes, and n to reflect partition after
        // sampling
        //
        // If record j is the only record in its new cluster from file k,
        // update inclusion patterns
        if(clust_sizes(k - 1, samp - 1) == 0){
            // Find the inclusion pattern of record j's cluster w/o record j
            arma::uvec ind = clust_sizes.col(samp - 1) > 0;
            // If j isn't in a singleton cluster, decrease the number of
            // clusters with the inclusion pattern of the cluster without
            // record j
            if(sum(ind) > 0){
                int oldind = bin_to_int_rcpp(ind, powers);
                cont(oldind)--;
            }
            // Update n if j is in a singleton cluster
            else{ n++; }
            // Increase the number of clusters with the inclusion pattern of
            // the cluster including record j
            ind(k - 1) = 1;
            int newind = bin_to_int_rcpp(ind, powers);
            cont(newind)++;
        }
        // Update clust_sizes
        clust_sizes(k - 1, samp - 1)++;
        clust_sizes_collapsed(samp - 1)++;
        // Using Z_members slows down inference when we have duplicates
        if(no_dups){
            // Update Z_members
            Z_members.submat(samp - 1, clust_sizes_collapsed(samp - 1) - 1,
                             samp - 1, clust_sizes_collapsed(samp - 1) - 1)
            = j + 1;
        }
    }

    return(List::create(Named("Z") = Z,
                        Named("cont") = cont,
                        Named("clust_sizes") = clust_sizes,
                        Named("n") = n,
                        Named("Z_members") = Z_members,
                        Named("clust_sizes_collapsed") =
                            clust_sizes_collapsed));
}

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
                              int num_field, const arma::vec& comparison_rps_probs){

    if(indexing_used){
        stop("Chaperones not implemented for use with indexing yet.");
    }

    // Perform num_chap_iter Chaperones updates
    for(int t = 0; t < num_chap_iter; t++){
        // Sample Chaperones
        int chaperone_1 = -99;
        int chaperone_2 = -99;
        // Uniform Chaperones distribution
        if(chap_type == 0){
            IntegerVector chaperones =
                uniform_chaperones(file_size_cum, valid_fp_matrix, fp_probs);
            chaperone_1 = min(chaperones) - 1;
            chaperone_2 = max(chaperones) - 1;
        }
        // Betancourt 2022 Chaperones distribution
        else if(chap_type == 1){
            IntegerVector chaperones =
                nonuniform_chaperones(file_size_cum, valid_fp_matrix, fp_probs,
                                      record_pairs, comparison_rps,
                                      comparison_rps_length, comparison_rps_probs);
            chaperone_1 = min(chaperones) - 1;
            chaperone_2 = max(chaperones) - 1;
        }

        // Get the label of the clusters the chaperones are in
        int label_1 = Z(chaperone_1);
        int label_2 = Z(chaperone_2);

        arma::uvec records;
        if(label_1 == label_2){
            records = arma::sort(arma::unique(arma::find(Z == label_1)));
        }
        else{
            records = arma::sort(arma::unique(arma::find(Z == label_1 || Z == label_2)));
        }

        for(int restrict_ind = 0; restrict_ind < num_restrict; restrict_ind++){
            for(int ind = 0; ind < records.size(); ind++){
                // Get the current record
                int j = records(ind);
                // Get the label of the cluster record j is currently in
                int label = Z(j);
                // k: the file number of the current record
                int k = file_labels(j);


                arma::vec valid_c = arma::vec(2);
                int no_move_flag = 0;

                if(j == chaperone_1 | j == chaperone_2){
                    if(label_1 == label_2){
                        // Two possible clusters, stay or make singleton
                        valid_c(0) = label;
                        valid_c(1) = -1;
                    }
                    else if(clust_sizes_collapsed(label - 1) == 1){
                        // Two possible clusters, stay (i.e. empty cluster) or join other cluster
                        // Need to check that joining other cluster is legal
                        if(j == chaperone_1){
                            if(clust_sizes(k - 1, label_2 - 1) < dup_upper_bound(k - 1)){
                                // Joining other cluster is legal
                                valid_c(0) = -1;
                                valid_c(1) = label_2;
                            }
                            else{
                                // No possible moves
                                no_move_flag = 1;
                            }
                        }
                        else if(j == chaperone_2){
                            if(clust_sizes(k - 1, label_1 - 1) < dup_upper_bound(k - 1)){
                                // Joining other cluster is legal
                                valid_c(0) = -1;
                                valid_c(1) = label_1;
                            }
                            else{
                                // No possible moves
                                no_move_flag = 1;
                            }
                        }
                        else{
                            Rcout << "Something went wrong, can't be here 1." << j << std::endl;
                        }
                    }
                    else if(clust_sizes_collapsed(label - 1) > 1){
                        // No possible moves
                        no_move_flag = 1;
                    }
                    else{
                        Rcout << "Something went wrong, can't be here 2." << j << std::endl;
                    }

                }
                else{
                    // Record is not a chaperone
                    if(label_1 == label_2){
                        // No possible moves
                        no_move_flag = 1;
                    }
                    else{
                        // Two possible clusters, stay or join other cluster
                        // Need to check that joining other cluster is legal
                        if(label == label_1){
                            if(clust_sizes(k - 1, label_2 - 1) < dup_upper_bound(k - 1)){
                                // Joining other cluster is legal
                                valid_c(0) = label_1;
                                valid_c(1) = label_2;
                            }
                            else{
                                // No possible moves
                                no_move_flag = 1;
                            }
                        }
                        else if(label == label_2){
                            if(clust_sizes(k - 1, label_1 - 1) < dup_upper_bound(k - 1)){
                                // Joining other cluster is legal
                                valid_c(0) = label_1;
                                valid_c(1) = label_2;
                            }
                            else{
                                // No possible moves
                                no_move_flag = 1;
                            }
                        }
                        else{
                            Rcout << "Something went wrong, can't be here 2." << j << std::endl;
                        }
                    }
                }

                if(no_move_flag == 0){
                    // Update cont, clust_sizes, and n to reflect partition without
                    // current record
                    //
                    // If record j is the only record in its cluster from file k, update
                    // inclusion patterns
                    if(clust_sizes(k - 1, label - 1) == 1){
                        // Find the inclusion pattern of record j's cluster
                        arma::uvec ind = clust_sizes.col(label - 1) > 0;
                        int oldind = bin_to_int_rcpp(ind, powers);

                        // Reduce the number of clusters with the inclusion pattern of
                        // the cluster including record j
                        cont(oldind)--;

                        // If j wasn't in a singleton cluster, increase the number of
                        // clusters with the inclusion pattern of the cluster without
                        // record j
                        if(sum(ind) > 1){
                            ind(k - 1) = 0;
                            int newind = bin_to_int_rcpp(ind, powers);
                            cont(newind)++;
                        }
                        // Update n if j was in a singleton cluster

                        else{ n--; }
                    }

                    // Using Z_members slows down inference when we have duplicates
                    if(no_dups){
                        int label_size = clust_sizes_collapsed(label - 1);
                        // If j was in a singleton cluster, update Z_members
                        if(label_size == 1){
                            Z_members.submat(label - 1, 0, label - 1, 0) = 0;
                        }
                        // If j wasn't in a singleton cluster, update Z_members
                        else if(label_size > 1){
                            arma::urowvec temp_1 = Z_members.submat(label - 1, 0,
                                                                    label - 1,
                                                                    label_size - 1);
                            temp_1.subvec(0, label_size - 2) =
                                temp_1.elem(arma::find(temp_1 != (j + 1))).t();
                            temp_1.subvec(label_size - 1, label_size - 1) = 0;
                            Z_members.submat(label - 1, 0,
                                             label - 1, label_size - 1) = temp_1;
                        }
                    }

                    // Update clust_sizes
                    clust_sizes(k - 1, label - 1)--;
                    clust_sizes_collapsed(label - 1)--;

                    // Give current record an empty label
                    Z(j) = -1;
                    arma::urowvec not_empty_clust;
                    if(no_dups){
                        not_empty_clust = clust_sizes_collapsed.t() > 0;
                    }

                    // Calculate cluster assignment probabilities
                    arma::vec log_probs = arma::zeros(valid_c.size());
                    for(int c_index = 0; c_index < valid_c.size(); c_index++){
                        // c_label: the label of the current cluster we're calculating
                        // the log_probs for
                        int c_label = valid_c(c_index);
                        // If c_label isn't an empty cluster
                        if(c_label != -1){
                            // Get the records in c_label
                            arma::uvec records;
                            // Using Z_members slows down inference when we have duplicates
                            if(no_dups){
                                records = (Z_members.submat(c_label - 1, 0, c_label - 1,
                                                            clust_sizes_collapsed(c_label
                                                                                      - 1)
                                                                - 1) - 1).t();
                            }
                            else{
                                records = arma::find(Z == c_label);
                            }

                            // Add the likelihood contribution
                            arma::uvec pairs = rp_ind.col(j);
                            pairs = pairs.elem(records) - 1;
                            if(sum(pairs < 0) > 0){
                                Rcout << "Something went wrong, not getting the correct\
                    record pairs" << j << std::endl;
                            }
                            log_probs(c_index) += arma::sum(log_like.elem(pairs));

                            if(!flat){
                                // c_k: the number of records from file k in the current
                                // cluster
                                int c_k = clust_sizes(k - 1, c_label - 1);
                                if(c_k >= dup_upper_bound(k - 1)){
                                    Rcout << "Something went wrong, trying to sample a\
                        cluster with too many records from file: " << k <<
                            std::endl;
                                }
                                // If there are no other records from file k in c_label
                                if(c_k == 0){
                                    arma::uvec ind = clust_sizes.col(c_label - 1) > 0;
                                    int ind_without = bin_to_int_rcpp(ind, powers);
                                    ind(k - 1) = 1;
                                    int ind_with = bin_to_int_rcpp(ind, powers);

                                    log_probs(c_index) +=
                                        as<arma::vec>(dup_count_prior[k - 1])(0) +
                                        std::log(cont(ind_with) + alphas(ind_with)) -
                                        std::log(cont(ind_without) + alphas(ind_without)
                                                     - 1);
                                }
                                // If there are other records from file k in c_label
                                else{
                                    log_probs(c_index) +=
                                        std::log(c_k + 1) +
                                        as<arma::vec>(dup_count_prior[k - 1])(c_k) -
                                        as<arma::vec>(dup_count_prior[k - 1])(c_k - 1);
                                }
                            }
                        }
                        // If c_label is an empty cluster
                        else{
                            if(!flat){
                                log_probs(c_index) +=
                                    as<arma::vec>(dup_count_prior(k - 1))(0) +
                                    std::log(n + 1) +
                                    std::log(cont(singleton_ind(k - 1) - 1) +
                                    alphas(singleton_ind(k - 1) - 1)) -
                                    std::log(n + alpha_0) + n_prior(n) - n_prior(n - 1);
                            }
                        }
                    }

                    //Sample the cluster assignments
                    arma::vec probs = exp(log_probs);

                    int samp = RcppArmadillo::sample(valid_c, 1, TRUE, probs)(0);
                    if(samp == -1){
                        if(no_dups){
                            arma::uvec empty_clusts =
                                arma::find((not_empty_clust == 0)) + 1;
                            samp = empty_clusts(0);
                        }
                        else{
                            for(int find = 1; find < r+1; find++){
                                arma::uvec ind = Z == find;
                                if(sum(ind) == 0){
                                    samp = find;
                                    find = r+1;
                                }
                            }
                        }
                    }
                    Z(j) = samp;
                    // Update cont, clust_sizes, and n to reflect partition after
                    // sampling
                    //
                    // If record j is the only record in its new cluster from file k,
                    // update inclusion patterns
                    if(clust_sizes(k - 1, samp - 1) == 0){
                        // Find the inclusion pattern of record j's cluster w/o record j
                        arma::uvec ind = clust_sizes.col(samp - 1) > 0;
                        // If j isn't in a singleton cluster, decrease the number of
                        // clusters with the inclusion pattern of the cluster without
                        // record j
                        if(sum(ind) > 0){
                            int oldind = bin_to_int_rcpp(ind, powers);
                            cont(oldind)--;
                        }
                        // Update n if j is in a singleton cluster
                        else{ n++; }
                        // Increase the number of clusters with the inclusion pattern of
                        // the cluster including record j
                        ind(k - 1) = 1;
                        int newind = bin_to_int_rcpp(ind, powers);
                        cont(newind)++;
                    }
                    // Update clust_sizes
                    clust_sizes(k - 1, samp - 1)++;
                    clust_sizes_collapsed(samp - 1)++;
                    // Using Z_members slows down inference when we have duplicates
                    if(no_dups){
                        // Update Z_members
                        Z_members.submat(samp - 1, clust_sizes_collapsed(samp - 1) - 1,
                                         samp - 1, clust_sizes_collapsed(samp - 1) - 1)
                        = j + 1;
                    }
                }

                if(j == chaperone_1){
                    label_1 = Z(j);
                }
                else if(j == chaperone_2){
                    label_2 = Z(j);
                }
            }
        }


    }
    if(extra_gibbs == 0){
        return(List::create(Named("Z") = Z,
                            Named("cont") = cont,
                            Named("clust_sizes") = clust_sizes,
                            Named("n") = n,
                            Named("Z_members") = Z_members,
                            Named("clust_sizes_collapsed") =
                                clust_sizes_collapsed));
    }
    else{
        return(sample_Z_rcpp(Z, clust_sizes, n, cont,
                             log_like, alphas,
                             alpha_0, dup_upper_bound,
                             dup_count_prior, n_prior, r,
                             r_1, valid_rp,
                             singleton_ind, rp_ind,
                             file_labels, powers,
                             flat, no_dups, cc, Z_members,
                             clust_sizes_collapsed, indexing_used));
    }

}
