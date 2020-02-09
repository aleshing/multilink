#include "sample_Z.h"
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

int bin_to_int_rcpp(const arma::uvec& bin, const arma::vec& powers){
    arma::vec temp = bin.t()*powers;
    return(temp(0));
}

List sample_Z_rcpp(arma::vec Z, arma::mat clust_sizes, int n, arma::vec cont,
                   const arma::vec& log_like, const arma::vec& alphas,
                   int alpha_0, const arma::vec& dup_upper_bound,
                   List dup_count_prior, const arma::vec& n_prior, int r,
                   int r_1, const arma::mat& valid_rp,
                   const arma::vec& singleton_ind, const arma::umat& rp_ind,
                   const arma::vec& file_labels, const arma::vec& powers,
                   int flat, int no_dups, int cc, arma::umat Z_members,
                   arma::vec clust_sizes_collapsed){
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
        clust_sizes(k-1, label-1)--;
        clust_sizes_collapsed(label-1)--;

        // Give current record an empty label
        Z(j) = -1;
        arma::vec valid_c;
        arma::urowvec not_empty_clust;
        if(no_dups){
            // Valid clusters are just clusters where there are already records,
            // where none of the records are from file k, along with the empty
            // cluster
            not_empty_clust = clust_sizes_collapsed.t()>0;
            arma::urowvec ind = clust_sizes.row(k-1)==0 && not_empty_clust;
            valid_c = arma::conv_to<arma::vec>::from(arma::find(ind.t())+1);
            int len = valid_c.size();
            valid_c.resize(len+1);
            valid_c(len) = -1;
        }
        else{
            // Find valid clusters, where -1 from the empty label will correspond
            // to forming a new cluster. Make sure to sort them!
            arma::vec possible_c = arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)))));
            arma::vec not_possible_c;
            // Have to do this next line since some of the clusters w/ valid records
            // may also contain non-valid records if we don't have transitive closures
            // of all connected components. Make sure to sort them!
            if(!cc){
                not_possible_c = arma::sort(arma::unique(Z.elem(arma::find(valid_rp.row(j)==0))));
            }

            // Can't join clusters that have more than the allowed number of
            // records from the same file! And make sure all of the clusters
            // in possible_c are actually possible if we don't have transitive
            // closures of all connected components
            if(!cc){
                int npc_counter = 0;
                int npc_size = not_possible_c.size();
                for(int i=0; i<possible_c.size(); i++){
                    int temp = possible_c(i);
                    if(temp!=-1){
                        if(clust_sizes(k-1, temp-1)>=dup_upper_bound(k-1)){
                            possible_c(i) = -1;
                        }
                        else if(npc_counter<npc_size){
                            int temp_npc=not_possible_c(npc_counter);
                            while(npc_counter<npc_size & temp_npc<=temp){
                                if(temp==temp_npc){
                                    possible_c(i) = -1;
                                }
                                npc_counter++;
                                if(npc_counter<npc_size){temp_npc=not_possible_c(npc_counter);}
                            }
                        }
                    }
                }
            }
            else{
                for(int i=0; i<possible_c.size(); i++){
                    int temp = possible_c(i);
                    if(temp!=-1){
                        if(clust_sizes(k-1, temp-1)>=dup_upper_bound(k-1)){
                            possible_c(i) = -1;
                        }
                    }
                }
            }
            valid_c = arma::unique(possible_c);
        }


        // Calculate cluster assignment probabilities
        arma::vec log_probs = arma::zeros(valid_c.size());
        for(int c_index=0; c_index<valid_c.size(); c_index++){
            // c_label: the label of the current cluster we're calculating
            // the log_probs for
            int c_label = valid_c(c_index);
            // If c_label isn't an empty cluster
            if(c_label!=-1){

                // Get the records in c_label
                arma::uvec records;
                // Using Z_members slows down inference when we have duplicates
                if(no_dups){
                    records = (Z_members.submat(c_label-1, 0, c_label-1,
                                                clust_sizes_collapsed(c_label-1)-1)-1).t();
                }
                else{
                    records = arma::find(Z==c_label);
                }

                // Add the likelihood contribution
                arma::uvec pairs = rp_ind.col(j);
                pairs = pairs.elem(records) - 1;
                if(sum(pairs<0)>0){
                    Rcout << "Something went wrong, not getting the correct\
                    record pairs" << j << std::endl;
                }
                log_probs(c_index) += arma::sum(log_like.elem(pairs));

                if(!flat){
                    // c_k: the number of records from file k in the current
                    // cluster
                    int c_k = clust_sizes(k-1, c_label-1);
                    if(c_k>=dup_upper_bound(k-1)){
                        Rcout << "Something went wrong, trying to sample a\
                        cluster with too many records from file: " << k << std::endl;
                    }
                    // If there are no other records from file k in c_label
                    if(c_k==0){
                        arma::uvec ind = clust_sizes.col(c_label-1)>0;
                        int ind_without = bin_to_int_rcpp(ind, powers);
                        ind(k-1) = 1;
                        int ind_with = bin_to_int_rcpp(ind, powers);

                        log_probs(c_index) +=
                            as<arma::vec>(dup_count_prior[k-1])(0) +
                            std::log(cont(ind_with) + alphas(ind_with)) -
                            std::log(cont(ind_without) + alphas(ind_without) - 1);
                    }
                    // If there are other records from file k in c_label
                    else{
                        log_probs(c_index) +=
                            std::log(c_k + 1) +
                            as<arma::vec>(dup_count_prior[k-1])(c_k) -
                            as<arma::vec>(dup_count_prior[k-1])(c_k-1);
                    }
                }
            }
            // If c_label is an empty cluster
            else{
                if(!flat){
                    log_probs(c_index) +=
                        as<arma::vec>(dup_count_prior(k-1))(0) +
                        std::log(n + 1) +
                        std::log(cont(singleton_ind(k-1)-1) + alphas(singleton_ind(k-1)-1)) -
                        std::log(n + alpha_0) + n_prior(n) - n_prior(n-1);
                }
            }
        }

        //Sample the cluster assignments
        arma::vec probs = exp(log_probs);

        int samp = RcppArmadillo::sample(valid_c, 1, TRUE, probs)(0);
        if(samp==-1){
            if(no_dups){
                arma::uvec empty_clusts = arma::find((not_empty_clust==0))+1;
                samp = empty_clusts(0);
            }
            else{
                for(int find=1; find<r+1; find++){
                    arma::uvec ind = Z==find;
                    if(sum(ind)==0){
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
        if(clust_sizes(k-1, samp-1)==0){
            // Find the inclusion pattern of record j's cluster w/o record j
            arma::uvec ind = clust_sizes.col(samp-1)>0;
            // If j isn't in a singleton cluster, decrease the number of
            // clusters with the inclusion pattern of the cluster without
            // record j
            if(sum(ind)>0){
                int oldind = bin_to_int_rcpp(ind, powers);
                cont(oldind)--;
            }
            // Update n if j is in a singleton cluster
            else{ n++; }
            // Increase the number of clusters with the inclusion pattern of
            // the cluster including record j
            ind(k-1) = 1;
            int newind = bin_to_int_rcpp(ind, powers);
            cont(newind)++;
        }
        // Update clust_sizes
        clust_sizes(k-1, samp-1)++;
        clust_sizes_collapsed(samp-1)++;

        // Using Z_members slows down inference when we have duplicates
        if(no_dups){
            // Update Z_members
            Z_members.submat(samp-1, clust_sizes_collapsed(samp-1)-1,
                             samp-1, clust_sizes_collapsed(samp-1)-1) = j+1;
        }
    }

    return(List::create(Named("Z")=Z,
                        Named("cont")=cont,
                        Named("clust_sizes")=clust_sizes,
                        Named("n")=n,
                        Named("Z_members")=Z_members,
                        Named("clust_sizes_collapsed")=clust_sizes_collapsed));
}
