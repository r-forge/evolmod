/*
 * tree_sum.h
 *
 *  Created on: Oct 7, 2013
 *      Author: hoblitz
 */

#ifndef TREE_SUM_H_
#define TREE_SUM_H_

#include "pclt.h"
#include "two_state_gibbs.h"
using namespace Rcpp;

// in-place functors for arma::transform()
struct ip_exp {
    ip_exp() {};
    double operator()(const double x) {
        return exp(x);
    };
};

struct ip_log {
    ip_log() {};
    double operator()(const double x) {
        return log(x);
    };
};

RcppExport SEXP treeConvolve(SEXP r_edge_matrix, SEXP r_branch_lengths,
        SEXP r_tip_states, SEXP r_rate_0, SEXP r_rate_1, SEXP r_n_max, SEXP r_root_p_0);

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

// annoying, cannot get arma::max to work with subviews
double colviewMax(const arma::subview_col<double> & col_view)
{
    const int len = col_view.n_elem;
    if (len == 0)
        Rcpp::stop("cannot take max of zero-length view");
    double m = col_view[0];
    for(int i = 1; i < len; i++) {
        if (m < col_view[i])
            m = col_view[i];
    }
    return m;
}

void updateRescaleValues(arma::cube & probs,
        double & rescale, const int & par_node_idx);

class TreeWorkspace {
private:
    arma::vec branch_lengths;
    arma::ivec tip_states;
    arma::imat edge_matrix;
    arma::cube q_probs, prior_probs, posterior_probs;
    double prior_rescale, post_rescale;
    int n_max, num_edges, num_tips;

    void initializeTips();

    inline void fillBranchProbs(const double & rate_01, const double & rate_10)
    {
        fillTransProb(0, rate_01, rate_10);
        fillTransProb(1, rate_10, rate_01);
    }

    void fillTransProb(const int & init_idx, const double & rate_init,
            const double & rate_next);

    void fillNodeProbs(const int & first_br_idx);

    double convolveBelowNode_0(const arma::cube & node_probs,
            const int & left_br_idx, const int & right_br_idx,
            const int & n) const;

    double convolveBelowNode_1(const arma::cube & node_probs,
            const int & left_br_idx, const int & right_br_idx,
            const int & n) const;

public:
    const arma::vec & getBranchLengths() const {
        return branch_lengths;
    }

    int getParentNode(const int & branch_idx) const {
        return edge_matrix(branch_idx, 0);
    }

    int getChildNode(const int & branch_idx) const {
        return edge_matrix(branch_idx, 1);
    }

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

    const arma::imat & getEdgeMatrix() const {
        return edge_matrix;
    }

    arma::mat convolveTree(const double & rate_01,
            const double & rate_10, const double & root_node_prob_0);

    TreeWorkspace(SEXP r_edge_matrix, SEXP r_branch_lengths, SEXP r_tip_states,
                SEXP r_n_max);
};

#endif /* TREE_SUM_H_ */
