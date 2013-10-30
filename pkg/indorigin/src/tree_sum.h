/*
 * tree_sum.h
 *
 *  Created on: Oct 7, 2013
 *      Author: hoblitz
 */

#ifndef TREE_SUM_H_
#define TREE_SUM_H_

#include "pclt.h"
using namespace Rcpp;
RcppExport SEXP treeCltSum(SEXP r_edge_matrix, SEXP r_branch_lengths,
        SEXP r_tip_states, SEXP r_rate_1, SEXP r_rate_2, SEXP r_n_max);


template <typename T>
void printElements(const arma::Mat<T> & A, const std::string & header) {
    typename arma::Mat<T>::const_iterator i;
    Rcout << header << std::endl;
    Rcout << "object dimensions = (" << A.n_rows << ", " << A.n_cols << ")";
    Rcout << std::endl;

    for(i = A.begin(); i != A.end(); i++) {
        Rcout << *i;
        Rcout << std::endl;
    }
    Rcout << std::endl;
    return;
}

template <typename T>
void printElements(const arma::Col<T> & A, const std::string & header) {
    typename arma::Col<T>::const_iterator i;
    Rcout << header << std::endl;
    Rcout << "object dimensions = (" << A.n_rows << ", " << A.n_cols << ")";
    Rcout << std::endl;

    for(i = A.begin(); i != A.end(); i++) {
        Rcout << *i;
        Rcout << std::endl;
    }
    Rcout << std::endl;
    return;
}

class TreeWorkspace {
private:
    arma::vec branch_lengths;
    arma::ivec tip_states;
    arma::imat edge_matrix;
    arma::cube q_probs, prior_probs, posterior_probs;
    arma::mat branch_frac_lik;

    int n_max, num_edges, num_tips, num_internal_nodes;
    double prior_rescale_value;
    double post_rescale_value;
    void initializeTips();
public:
    double getBranchLength(const int & idx) const {
        return branch_lengths[idx];
    }

    int getNumEdges() const {
        return num_edges;
    }

    const arma::cube & getQProbs() const{
        return q_probs;
    }

    const arma::cube & getPostProbs() const {
        return posterior_probs;
    }

    const arma::cube & getPriorProbs() const {
        return prior_probs;
    }

    TreeWorkspace(SEXP r_edge_matrix, SEXP r_branch_lengths, SEXP r_tip_states,
                SEXP r_n_max);

    void fillTransProb(const int & init_idx, const double & rate_init,
            const double & rate_next);

    void fillNodeProbs(const int & first_br_idx);

    void computeRootDist();

    double convolveBelowNode_0(const arma::cube & node_probs,
            const int & left_br_idx, const int & right_br_idx,
            const int & n) const;

    double convolveBelowNode_1(const arma::cube & node_probs,
            const int & left_br_idx, const int & right_br_idx,
            const int & n) const;
};

inline void fillBranchProbs(TreeWorkspace & tree, const double & rate_0, const double & rate_1);

#endif /* TREE_SUM_H_ */
