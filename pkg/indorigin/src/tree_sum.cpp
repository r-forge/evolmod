/*
 * tree_sum.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: hoblitz
 */

#include "tree_sum.h"

SEXP treeConvolve(SEXP r_edge_matrix, SEXP r_branch_lengths,
        SEXP r_tip_states, SEXP r_rate_0, SEXP r_rate_1, SEXP r_n_max, SEXP r_root_p_0)
{
    TreeWorkspace tree(r_edge_matrix, r_branch_lengths,
            r_tip_states, r_n_max);
    const double rate_0 = Rcpp::as<double>(r_rate_0);
    const double rate_1 = Rcpp::as<double>(r_rate_1);
    const double root_p_0 = Rcpp::as<double>(r_root_p_0);

    arma::mat output = tree.convolveTree(rate_0, rate_1, root_p_0);
    return wrap(output);
}

TreeWorkspace::TreeWorkspace(SEXP r_edge_matrix, SEXP r_branch_lengths, SEXP r_tip_states,
        SEXP r_n_max) : edge_matrix(Rcpp::as<arma::imat>(r_edge_matrix)),
                branch_lengths(Rcpp::as<arma::vec>(r_branch_lengths)),
                tip_states(Rcpp::as<arma::ivec>(r_tip_states)),
                post_rescale(0.0), prior_rescale(0.0),
                n_max(Rcpp::as<int>(r_n_max))
{
    edge_matrix -= 1;   //convert edge indices to C-style
    const int num_internal_nodes = edge_matrix.n_rows / 2;
    num_tips = tip_states.n_elem;
    num_edges = edge_matrix.n_rows;

    q_probs = arma::cube(2 * (n_max + 1) + 2, num_edges, 2, arma::fill::zeros);
    prior_probs = arma::cube(n_max + 1, num_tips + num_internal_nodes, 2, arma::fill::zeros);
    posterior_probs = arma::cube(n_max + 1, num_tips + num_internal_nodes, 2, arma::fill::zeros);
    initializeTips();
}

void TreeWorkspace::fillTransProb(const int & init_idx,
        const double & rate_init, const double & rate_next) {
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
            slice_view.col(edge_idx).transform(ip_exp());
        }
        prev_branch_len = cur_branch_len;
    }
    return;
}

inline void TreeWorkspace::initializeTips() {
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

arma::mat TreeWorkspace::convolveTree(const double & rate_01,
        const double & rate_10, const double & root_node_prob_0)
{
    fillBranchProbs(rate_01, rate_10);
    for(int i = 0; i < num_edges; i += 2) {
        fillNodeProbs(i);
    }
    // root node is now a rescaled probability
    const int root_node_idx = num_tips;
    IntegerVector tip_states_copy(tip_states.n_elem);
    for(int i = 0; i < tip_states.n_elem; i++) {
        tip_states_copy[i] = tip_states[i];
    }

    arma::vec arma_root_dist(2);
    arma_root_dist(0) = root_node_prob_0;
    arma_root_dist(1) = 1.0 - root_node_prob_0;

    // get likelihood of tip data
    const double post_lik = TwoStatePhyloLikelihood(edge_matrix + 1, tip_states_copy,
            branch_lengths, rate_01, rate_10, arma_root_dist);
    post_rescale -= log(post_lik);

    // gather root node vectors
    arma::mat output(n_max + 1, 2);
    arma::mat posterior_output =
            arma::join_rows(posterior_probs.slice(0).col(root_node_idx),
                    posterior_probs.slice(1).col(root_node_idx));
    arma::mat prior_output =
            arma::join_rows(prior_probs.slice(0).col(root_node_idx),
                    prior_probs.slice(1).col(root_node_idx));

    // transform output to log scale, rescale, and combine
    posterior_output.transform(ip_log());
    prior_output.transform(ip_log());

    posterior_output.col(0) += post_rescale + log(root_node_prob_0);
    posterior_output.col(1) += post_rescale + log(1.0 - root_node_prob_0);

    prior_output.col(0) += prior_rescale + log(root_node_prob_0);
    prior_output.col(1) += prior_rescale + log(1.0 - root_node_prob_0);

    for(int i = 0; i < output.n_rows; i++){
        output(i, 0) = logspaceAdd(prior_output(i, 0), prior_output(i, 1));
        output(i, 1) = logspaceAdd(posterior_output(i, 0), posterior_output(i, 1));
    }

    return output;
}

void TreeWorkspace::fillNodeProbs(const int & post_order_idx)
{
    const int left_br_idx = post_order_idx;
    const int right_br_idx = post_order_idx + 1;
    const int par_node_idx = getParentNode(post_order_idx);

    for(int n = 0; n <= n_max; n++) {
        posterior_probs(n, par_node_idx, 0) =
                convolveBelowNode_0(posterior_probs, left_br_idx, right_br_idx, n);

        posterior_probs(n, par_node_idx, 1) =
                convolveBelowNode_1(posterior_probs, left_br_idx, right_br_idx, n);

        prior_probs(n, par_node_idx, 0) =
                convolveBelowNode_0(prior_probs, left_br_idx, right_br_idx, n);

        prior_probs(n, par_node_idx, 1) =
                convolveBelowNode_1(prior_probs, left_br_idx, right_br_idx, n);
    }

    updateRescaleValues(prior_probs, prior_rescale, par_node_idx);
    updateRescaleValues(posterior_probs, post_rescale, par_node_idx);

    return;
}

void updateRescaleValues(arma::cube & probs,
        double & rescale, const int & par_node_idx)
{
    const double m0 = colviewMax(probs.slice(0).col(par_node_idx));
    const double m1 = colviewMax(probs.slice(1).col(par_node_idx));
    const double m = fmax(m0, m1);
    if (m <= 0.0)
        return;
    probs.tube(arma::span(), arma::span(par_node_idx)) /= m;
    rescale += log(m);
}

/* convolveBelowNode_(0/1) differ subtly according to how
 * the q_prob member is traversed, because we are summing forward
 * transitions. Unfortunately we must write two very similar functions.
 * Changing the storage of q_prob could alleviate this.
 */

double TreeWorkspace::convolveBelowNode_0(const arma::cube & node_probs,
        const int & left_br_idx, const int & right_br_idx,
        const int & n) const
{
    const int left_ch_node = getChildNode(left_br_idx);
    const int right_ch_node = getChildNode(right_br_idx);
    const int par_node = getParentNode(left_br_idx);
    if (par_node != getParentNode(right_br_idx)) {
        Rcpp::stop("left/right branches must have same parent");
    }

    double sum_left(0.0);
    double sum_right(0.0);
    double sum_total(0.0);
    int j, k;

    for(int i = 0; i <= n; i++) {
        /* sum over combinations of left node terminating
         * in state 0, starting in state 0 */
        for(j = 0, k = i; j <= i; j++, k--) {
            sum_left += node_probs(j, left_ch_node, 0) *
                    q_probs(2 * k, left_br_idx, 0);
        }

        /* sum over combinations of left node terminating in
         * state 1, starting in state 0 */
        for(j = i - 1, k = 1; j >= 0; j--, k++) {
            sum_left += node_probs(j, left_ch_node, 1) *
                    q_probs(2 * k - 1, left_br_idx, 0);
        }

        /* right node terminating in state 0, starting in state 0 */
        for(j = 0, k = n - i; j <= n - i; j++, k--) {
            sum_right += node_probs(j, right_ch_node, 0) *
                    q_probs(2 * k, right_br_idx, 0);
        }

        /* right node terminating in state 1, starting in state 0 */
        for(j = n - i - 1, k = 1; j >= 0; j--, k++) {
            sum_right += node_probs(j, right_ch_node, 1) *
                    q_probs(2 * k - 1, right_br_idx, 0);
        }
        sum_total += sum_left * sum_right;
        sum_left = 0.0;
        sum_right = 0.0;
    }
    return sum_total;
}

double TreeWorkspace::convolveBelowNode_1(const arma::cube & node_probs,
        const int & left_br_idx, const int & right_br_idx,
        const int & n) const
{
    const int left_ch_node = getChildNode(left_br_idx);
    const int right_ch_node = getChildNode(right_br_idx);
    const int par_node = getParentNode(left_br_idx);
    if (par_node != getParentNode(right_br_idx)) {
        Rcpp::stop("left/right branches must have same parent");
    }

    double sum_left(0.0);
    double sum_right(0.0);
    double sum_total(0.0);
    int j, k;

    for(int i = 0; i <= n; i++) {
        /* right node starting in state 1, ending in state 1 */
        for(j = 0, k = i; j <= i; j++, k--) {
            sum_left += node_probs(j, left_ch_node, 1) *
                    q_probs(2 * k, left_br_idx, 1);
        }
        /* left node starting in state 1, ending in state 0 */
        for(j = 0, k = i; j <= i; j++, k--) {
            sum_left += node_probs(j, left_ch_node, 0) *
                    q_probs(2 * k + 1, left_br_idx, 1);
        }
        /* right node starting in state 1, ending in state 0 */
        for(j = 0, k = n - i; j <= n - i; j++, k--) {
            sum_right += node_probs(j, right_ch_node, 0) *
                    q_probs(2 * k + 1, right_br_idx, 1);
        }
        /* right node starting in state 1, ending in state 1 */
        for(j = 0, k = n - i; j <= n - i; j++, k--) {
            sum_right += node_probs(j, right_ch_node, 1) *
                    q_probs(2 * k, right_br_idx, 1);
        }

        sum_total += sum_left * sum_right;
        sum_left = 0.0;
        sum_right = 0.0;
    }
    return sum_total;
}

