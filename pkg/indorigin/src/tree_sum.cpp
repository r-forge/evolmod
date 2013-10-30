/*
 * tree_sum.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: hoblitz
 */

#include "tree_sum.h"

SEXP treeCltSum(SEXP r_edge_matrix, SEXP r_branch_lengths,
        SEXP r_tip_states, SEXP r_rate_0, SEXP r_rate_1, SEXP r_n_max) {
    SEXP ret = R_NilValue;
    TreeWorkspace tree(r_edge_matrix, r_branch_lengths,
            r_tip_states, r_n_max);
    double rate_0 = Rcpp::as<double>(r_rate_0);
    double rate_1 = Rcpp::as<double>(r_rate_1);
    fillBranchProbs(tree, rate_0, rate_1);
    tree.getPostProbs().print();
    tree.getPriorProbs().print();
    tree.computeRootDist();
    return Rcpp::List::create(
            Named("prior") = tree.getPriorProbs(),
            Named("posterior") = tree.getPostProbs());
}

struct ip_exp { // in-place exp functor for arma::transform()
    ip_exp() {};
    double operator()(const double x) {
        return exp(x);
    };
};

TreeWorkspace::TreeWorkspace(SEXP r_edge_matrix, SEXP r_branch_lengths, SEXP r_tip_states,
        SEXP r_n_max) {
    edge_matrix = Rcpp::as<arma::imat>(r_edge_matrix);
    edge_matrix -= 1;

    branch_lengths = Rcpp::as<arma::vec>(r_branch_lengths);
    tip_states = Rcpp::as<arma::ivec>(r_tip_states);

    post_rescale_value = 0.0;
    prior_rescale_value = 0.0;
    n_max = Rcpp::as<int>(r_n_max);
    num_internal_nodes = edge_matrix.n_rows / 2;
    num_tips = tip_states.n_elem;
    num_edges = edge_matrix.n_rows;

    q_probs = arma::cube(2 * (n_max + 1) + 2, num_edges, 2, arma::fill::zeros);
    branch_frac_lik = arma::mat(n_max + 1, 2, arma::fill::zeros);
    prior_probs = arma::cube(n_max + 1, num_tips + num_internal_nodes, 2, arma::fill::zeros);
    posterior_probs = arma::cube(n_max + 1, num_tips + num_internal_nodes, 2, arma::fill::zeros);
    initializeTips();
}

void TreeWorkspace::fillTransProb(const int & init_idx,
        const double & rate_init, const double & rate_next) {
    ip_exp exp_functor;
    double cur_branch_len, prev_branch_len = -1.0;
    arma::mat & slice_view = q_probs.slice(init_idx);
    for(int edge_idx = 0; edge_idx < num_edges; edge_idx++) {
        cur_branch_len = getBranchLength(edge_idx);
        if (cur_branch_len == prev_branch_len) {
            // Branch lengths are same. No need to re-compute, just copy.
            slice_view.col(edge_idx) = slice_view.col(edge_idx - 1);
        }
        else {
            pclt(cur_branch_len, rate_init, rate_next, slice_view.col(edge_idx));
            slice_view.col(edge_idx).transform(exp_functor);
        }
        prev_branch_len = cur_branch_len;
    }
    return;
}

void TreeWorkspace::initializeTips() {
    prior_probs(arma::span(0), arma::span(0, num_tips - 1), arma::span()).fill(1.0);

    arma::cube::slice_iterator j0 = posterior_probs.begin_slice(0);
    arma::cube::slice_iterator j1 = posterior_probs.begin_slice(1);
    for(arma::ivec::const_iterator i = tip_states.begin();
            i != tip_states.end(); i++, j0+= (n_max + 1), j1+= (n_max + 1)) {
        if (*i == 0)
            *j0 = 1.0;
        else
            *j1 = 1.0;
    }
}

inline void fillBranchProbs(TreeWorkspace & tree, const double & rate_0, const double & rate_1) {
    tree.fillTransProb(0, rate_0, rate_1);
    tree.fillTransProb(1, rate_1, rate_0);
}

void TreeWorkspace::computeRootDist() {
    for(int i = 0; i < num_edges; i += 2) {
        fillNodeProbs(i);
    }
}

void TreeWorkspace::fillNodeProbs(const int & post_order_idx) {
    const int left_br_idx = post_order_idx;
    const int right_br_idx = post_order_idx + 1;
    const int par_node_idx = edge_matrix(post_order_idx, 0);
    double prob;
    Rcout << "fill node probs: " << post_order_idx <<
            " left: " << left_br_idx <<
            " right: " << right_br_idx << std::endl;
    for(int i = 0; i <= n_max; i++) {
        Rcout << "    i = " << i << std::endl;
        prob = convolveBelowNode_0(posterior_probs,
                left_br_idx, right_br_idx, i);
        posterior_probs(i, par_node_idx, 0) = prob;

        Rcout << "        node 0 done ;" << std::endl;

        prob = convolveBelowNode_1(posterior_probs,
                left_br_idx, right_br_idx, i);
        posterior_probs(i, par_node_idx, 1) = prob;

        Rcout << "        node 1 done ;" << std::endl;

        prob = convolveBelowNode_0(prior_probs,
                left_br_idx, right_br_idx, i);

        prior_probs(i, par_node_idx, 0) = prob;
        prob = convolveBelowNode_1(prior_probs,
                left_br_idx, right_br_idx, i);
        prior_probs(i, par_node_idx, 1) = prob;
    }

    return;
}

/* convolveBelowNode_(0/1) differ only slightly according to how
 * the q_prob member is traversed, because we are summing forward
 * transitions. Unfortunately we must write two very similar functions.
 */

double TreeWorkspace::convolveBelowNode_0(const arma::cube & node_probs,
        const int & left_br_idx, const int & right_br_idx,
        const int & n) const
{
    const int par_node = edge_matrix(left_br_idx, 0);
    if (par_node != edge_matrix(right_br_idx, 0)) {
        Rcpp::stop("left/right branches must have same parent");
    }
    const int left_ch_node = edge_matrix(left_br_idx, 1);
    const int right_ch_node = edge_matrix(right_br_idx, 1);
    double acc_left(0.0);
    double acc_right(0.0);
    double acc_total(0.0);

    for(int i = 0; i <= n; i++) {
        for(int j = 0, k = i; j <= i; j++, k--) {
            acc_left += node_probs(j, left_ch_node, 0) *
                    q_probs(2 * k, left_br_idx, 0);
        }

        for(int j = i - 1, k = 1; j >= 0; j--, k++) {
            acc_left += node_probs(j, left_ch_node, 1) *
                    q_probs(2 * k - 1, left_br_idx, 0);
        }

        for(int j = 0, k = n - i; j <= n - i; j++, k--) {
            acc_right += node_probs(j, right_ch_node, 0) *
                    q_probs(2 * k, right_br_idx, 0);
        }

        for(int j = n - i - 1, k = 1; j >= 0; j--, k++) {
            acc_right += node_probs(j, right_ch_node, 1) *
                    q_probs(2 * k - 1, right_br_idx, 0);
        }
        acc_total += acc_left * acc_right;
        acc_left = 0.0;
        acc_right = 0.0;
    }
    return acc_total;
}

double TreeWorkspace::convolveBelowNode_1(const arma::cube & node_probs,
        const int & left_br_idx, const int & right_br_idx,
        const int & n) const
{
    const int par_node = edge_matrix(left_br_idx, 0);
    if (par_node != edge_matrix(right_br_idx, 0)) {
        Rcpp::stop("left/right branches must have same parent");
    }
    const int left_ch_node = edge_matrix(left_br_idx, 1);
    const int right_ch_node = edge_matrix(right_br_idx, 1);
    double acc_left(0.0);
    double acc_right(0.0);
    double acc_total(0.0);

    for(int i = 0; i <= n; i++) {
        for(int j = 0, k = i; j <= i; j++, k--) {
            acc_left += node_probs(j, left_ch_node, 1) *
                    q_probs(2 * k, left_br_idx, 1);
        }

        for(int j = 0, k = i; j <= i; j++, k--) {
            acc_left += node_probs(j, left_ch_node, 0) *
                    q_probs(2 * k + 1, left_br_idx, 1);
        }

        for(int j = 0, k = n - i; j <= n - i; j++, k--) {
            acc_right += node_probs(j, right_ch_node, 0) *
                    q_probs(2 * k, right_br_idx, 1);
        }

        for(int j = 0, k = n - i; j <= n - i; j++, k--) {
            acc_right += node_probs(j, right_ch_node, 1) *
                    q_probs(2 * k + 1, right_br_idx, 1);
        }

        acc_total += acc_left * acc_right;
        acc_left = 0.0;
        acc_right = 0.0;
    }
    return acc_total;
}









































