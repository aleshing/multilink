#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec get_same_part(int TT, int r, const arma::mat& partitions){

    // Find the first partition in the chain that's the same as the current
    // partition (if it exists)
    arma::vec same_part = arma::zeros<arma::vec>(TT) - 1;
    for(int t = 1; t < TT; t++){
        for(int l = 0; l < t; l++){
            if(sum(partitions.col(t) == partitions.col(l)) == r){
                same_part(t) = l;
                // Short the loop
                l = t + 1;
            }
        }
    }

    return(same_part);
}

arma::mat get_delta_mats(int TT, int r, const arma::mat& partitions,
                         const arma::vec& same_part){

    arma::mat delta_mats = arma::zeros<arma::mat>(r * TT, r);
    // Calculate the similarity matrix for each sample
    for(int t = 0; t < TT; t++){
        checkUserInterrupt();
        if(same_part(t) > -1){
            delta_mats.submat(r * t, 0, r * (t + 1) - 1, r - 1) =
                delta_mats.submat(r * same_part(t), 0,
                                  r * (same_part(t) + 1) - 1, r - 1);
        }
        else{
            for(int i = 0; i < r; i++){
                for(int j = 0; j < i; j++){
                    if(partitions(i, t) == partitions(j, t)){
                        delta_mats(i + r * t, j) = 1;
                        delta_mats(j + r * t, i) = 1;
                    }
                }
                delta_mats(i + r * t, i) = 1;
            }
        }
    }
    return(delta_mats);
}

// Third version
arma::rowvec get_FM2_loss(int TT, int r, const arma::umat& running,
                            double L_FM2, const arma::mat& delta_mats,
                            const arma::mat& A, const arma::vec& same_part){
    arma::rowvec post_loss_FM2 = arma::zeros<arma::rowvec>(TT);

    for(int i = 0; i < r; i++){
        arma::mat temp_loss = arma::zeros<arma::mat>(TT, TT);
        arma::mat temp_prod = arma::zeros<arma::mat>(TT, TT);
        for(int t = 0; t < TT; t++){
            if(!running(i, t) & (A(i, t) != r - 1)){
                if(same_part(t) > -1){
                    temp_loss.row(t) = temp_loss.row(same_part(t));
                }
                else{
                    checkUserInterrupt();

                    arma::rowvec delta_temp = delta_mats.row(i + r * t);
                    for(int l = 0; l < TT; l++){
                        if(!running(i, l) & (t != l)){
                            if(same_part(l) > -1){
                                temp_loss(t, l) = temp_loss(t, same_part(l));
                            }
                            else if(l < t){
                                if((A(i, l) - temp_prod(l, t)) > -1){
                                    temp_loss(t, l) = L_FM2 / TT;
                                }
                            }
                            else{
                                double temp_temp_prod =
                                    dot(delta_mats.row(i + r * l), delta_temp);
                                temp_prod(t, l) = temp_temp_prod;
                                if((A(i, l) - temp_temp_prod) > -1){
                                    temp_loss(t, l) = L_FM2 / TT;
                                }
                            }
                        }
                    }
                }
            }
        }
        post_loss_FM2 += sum(temp_loss, 1).t();
    }
    return(post_loss_FM2);
}

arma::rowvec get_FM2_abstain(int TT, int r, const arma::umat& running,
                             double L_FM2, const arma::mat& delta_mats,
                             const arma::mat& A, const arma::vec& same_part,
                             const arma::vec& power_row,
                             const arma::mat& new_delta_mats,
                             const arma::umat& new_running,
                             const arma::mat& new_A){
    arma::rowvec post_loss_FM2 = arma::zeros<arma::rowvec>(TT);

    for(int i = 0; i < r; i++){
        // If record i isn't abstaining
        if(power_row(i) == 1){
            arma::mat temp_loss = arma::zeros<arma::mat>(TT, TT);
            for(int t = 0; t < TT; t++){
                if(!new_running(i, t)){
                    if(same_part(t) > -1){
                        temp_loss.row(t) = temp_loss.row(same_part(t));
                    }
                    else{
                        checkUserInterrupt();

                        arma::rowvec new_delta_temp =
                            new_delta_mats.row(i + r * t);
                        for(int l = 0; l < TT; l++){
                            if(!running(i, l) & (t != l)){
                                if(same_part(l) > -1){
                                    temp_loss(t, l) =
                                        temp_loss(t, same_part(l));
                                }
                                else if((A(i, l) -
                                        dot(delta_mats.row(i + r * l),
                                            new_delta_temp)) > -1){
                                    temp_loss(t, l)= L_FM2 / TT;
                                }
                            }
                        }
                    }
                }
            }
            post_loss_FM2 += sum(temp_loss, 1).t();
        }
    }
    return(post_loss_FM2);
}

arma::rowvec get_post_loss_abstain_row(int TT, int r, double L_FNM,
                                       double L_FM1, double L_FM2, double L_A,
                                       const arma::vec& post_prob_match,
                                       const arma::vec&post_prob_no_match,
                                       const arma::umat& running,
                                       const arma::vec& power_row,
                                       const arma::vec& same_part,
                                       const arma::mat& partitions,
                                       const arma::mat& delta_mats,
                                       const arma::mat& A){
    arma::mat new_delta_mats = delta_mats;
    for(int i = 0; i < r; i++){
        if(power_row(i) == 0){
            new_delta_mats.col(i) = arma::zeros<arma::vec>(TT * r);
            for(int t = 0; t < TT; t++){
                new_delta_mats.row(i + r * t) = arma::zeros<arma::rowvec>(r);
                new_delta_mats(i + r * t, i) = 1;
            }
        }
    }

    arma::umat new_running = arma::umat(r, TT);
    arma::mat new_A = arma::zeros<arma::mat>(r, TT);
    for(int t = 0; t < TT; t++){
        if(same_part(t) > -1){
            new_A.col(t) = new_A.col(same_part(t));
            new_running.col(t) = new_running.col(same_part(t));
        }
        else{
            new_A.col(t) = sum(new_delta_mats.rows(r * t, r * (t + 1) - 1), 1)
            - 1;
            new_running.col(t) = new_A.col(t) == 0;
        }
    }

    arma::rowvec new_post_loss = arma::zeros<arma::rowvec>(TT);

    new_post_loss += L_FNM * post_prob_match.t() *
        (arma::conv_to<arma::mat>::from(new_running).each_col() % power_row);
    new_post_loss += L_FM1 * post_prob_no_match.t() *
        ((arma::conv_to<arma::mat>::from(1 - new_running).each_col())
             % power_row);
    new_post_loss += L_A * sum(1 - power_row);

    new_post_loss += get_FM2_abstain(TT, r, running, L_FM2, delta_mats, A,
                                     same_part, power_row, new_delta_mats,
                                     new_running, new_A);
    return(new_post_loss);
}

// [[Rcpp::export]]
arma::rowvec get_posterior_loss_allcpp(int TT, int r,
                                       const arma::mat& partitions,
                                       double L_FNM, double L_FM1,
                                       double L_FM2){
    arma::vec same_part = get_same_part(TT, r, partitions);
    arma::mat delta_mats = get_delta_mats(TT, r, partitions, same_part);

    // A(i, t) = \sum_{j\neq i}\delta^t_{ij}
    // running(i, t) = I(\sum_{j\neq i}\delta^t_{ij}=0)
    arma::umat running = arma::umat(r, TT);
    arma::mat A = arma::zeros<arma::mat>(r, TT);
    for(int t = 0; t < TT; t++){
        if(same_part(t) > -1){
            A.col(t) = A.col(same_part(t));
            running.col(t) = running.col(same_part(t));
        }
        else{
            A.col(t) = sum(delta_mats.rows(r * t, r * (t + 1) - 1), 1) - 1;
            running.col(t) = A.col(t) == 0;
        }
    }

    // Precompute P(\sum_{j\neq i}\delta_{ij}=0|-)
    arma::colvec post_prob_no_match =
        sum(arma::conv_to<arma::mat>::from(running), 1) / TT;
    arma::vec post_prob_match = 1 - post_prob_no_match;

    arma::rowvec post_loss = arma::zeros<arma::rowvec>(TT);
    // Calculate the loss, for records in the partitions that made no matches
    post_loss += L_FNM * post_prob_match.t() * running;
    // Calculate part of the loss, for records in the partitions that made
    // matches
    post_loss += L_FM1 * post_prob_no_match.t() * (1 - running);

    post_loss += get_FM2_loss(TT, r, running, L_FM2, delta_mats, A, same_part);

    return(post_loss);
}

// [[Rcpp::export]]
arma::rowvec get_posterior_loss_abstain_cpp(int TT, int r,
                                            const arma::mat& partitions,
                                            double L_FNM, double L_FM1,
                                            double L_FM2, double L_A,
                                            const arma::mat& power){

    arma::vec same_part = get_same_part(TT, r, partitions);
    arma::mat delta_mats = get_delta_mats(TT, r, partitions, same_part);

    // running(i, t) = I(\sum_{j\neq i}\delta^t_{ij}=0)
    arma::umat running = arma::umat(r, TT);
    arma::mat A = arma::zeros<arma::mat>(r, TT);
    for(int t = 0; t < TT; t++){
        if(same_part(t) > -1){
            A.col(t) = A.col(same_part(t));
            running.col(t) = running.col(same_part(t));
        }
        else{
            A.col(t) = sum(delta_mats.rows(r * t, r * (t + 1) - 1), 1) - 1;
            running.col(t) = A.col(t) == 0;
        }
    }

    // Precompute P(\sum_{j\neq i}\delta_{ij}=0|-)
    arma::colvec post_prob_no_match =
        sum(arma::conv_to<arma::mat>::from(running), 1) / TT;
    arma::vec post_prob_match = 1 - post_prob_no_match;

    arma::rowvec post_loss = arma::zeros<arma::rowvec>(TT * power.n_rows);
    // Calculate the loss, for records in the partitions that made no matches
    post_loss.head(TT) += L_FNM * post_prob_match.t() * running;
    // Calculate part of the loss, for records in the partitions that made matches
    post_loss.head(TT) += L_FM1 * post_prob_no_match.t() * (1 - running);

    post_loss.head(TT) += get_FM2_loss(TT, r, running, L_FM2, delta_mats, A,
                   same_part);

    for(int j = 1; j < power.n_rows; j++){
        post_loss.subvec(TT * j, TT * (j + 1) - 1) =
            get_post_loss_abstain_row(TT, r, L_FNM, L_FM1, L_FM2, L_A,
                                      post_prob_match, post_prob_no_match,
                                      running, power.row(j).t(), same_part,
                                      partitions, delta_mats, A);
    }

    return(post_loss);
}

