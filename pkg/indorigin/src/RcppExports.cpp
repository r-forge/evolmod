// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// twoStateSufficientStatistics
NumericVector twoStateSufficientStatistics(arma::Mat<int>& treeEdges, IntegerVector& tipStates, NumericVector& branchLengths, double lambda_01, double lambda_10, NumericVector& rootDist);
RcppExport SEXP indorigin_twoStateSufficientStatistics(SEXP treeEdgesSEXP, SEXP tipStatesSEXP, SEXP branchLengthsSEXP, SEXP lambda_01SEXP, SEXP lambda_10SEXP, SEXP rootDistSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::Mat<int>& >::type treeEdges(treeEdgesSEXP );
        Rcpp::traits::input_parameter< IntegerVector& >::type tipStates(tipStatesSEXP );
        Rcpp::traits::input_parameter< NumericVector& >::type branchLengths(branchLengthsSEXP );
        Rcpp::traits::input_parameter< double >::type lambda_01(lambda_01SEXP );
        Rcpp::traits::input_parameter< double >::type lambda_10(lambda_10SEXP );
        Rcpp::traits::input_parameter< NumericVector& >::type rootDist(rootDistSEXP );
        NumericVector __result = twoStateSufficientStatistics(treeEdges, tipStates, branchLengths, lambda_01, lambda_10, rootDist);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// twoStateCompleteDataLogPosterior
double twoStateCompleteDataLogPosterior(NumericVector suffStat, double lambda_01, double lambda_10, double prior_alpha_01, double prior_beta_01, double prior_alpha_10, double prior_beta_10);
RcppExport SEXP indorigin_twoStateCompleteDataLogPosterior(SEXP suffStatSEXP, SEXP lambda_01SEXP, SEXP lambda_10SEXP, SEXP prior_alpha_01SEXP, SEXP prior_beta_01SEXP, SEXP prior_alpha_10SEXP, SEXP prior_beta_10SEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type suffStat(suffStatSEXP );
        Rcpp::traits::input_parameter< double >::type lambda_01(lambda_01SEXP );
        Rcpp::traits::input_parameter< double >::type lambda_10(lambda_10SEXP );
        Rcpp::traits::input_parameter< double >::type prior_alpha_01(prior_alpha_01SEXP );
        Rcpp::traits::input_parameter< double >::type prior_beta_01(prior_beta_01SEXP );
        Rcpp::traits::input_parameter< double >::type prior_alpha_10(prior_alpha_10SEXP );
        Rcpp::traits::input_parameter< double >::type prior_beta_10(prior_beta_10SEXP );
        double __result = twoStateCompleteDataLogPosterior(suffStat, lambda_01, lambda_10, prior_alpha_01, prior_beta_01, prior_alpha_10, prior_beta_10);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// twoStatePhyloGibbsSampler
NumericVector twoStatePhyloGibbsSampler(IntegerVector treeEdges, IntegerVector cubeDims, NumericMatrix branchLengths, NumericVector rootDist, IntegerMatrix tipStates, double initial_lambda_01, double initial_lambda_10, double prior_alpha_01, double prior_beta_01, double prior_alpha_10, double prior_beta_10, int mcmcSize, int mcmcBurnin, int mcmcSubsample);
RcppExport SEXP indorigin_twoStatePhyloGibbsSampler(SEXP treeEdgesSEXP, SEXP cubeDimsSEXP, SEXP branchLengthsSEXP, SEXP rootDistSEXP, SEXP tipStatesSEXP, SEXP initial_lambda_01SEXP, SEXP initial_lambda_10SEXP, SEXP prior_alpha_01SEXP, SEXP prior_beta_01SEXP, SEXP prior_alpha_10SEXP, SEXP prior_beta_10SEXP, SEXP mcmcSizeSEXP, SEXP mcmcBurninSEXP, SEXP mcmcSubsampleSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type treeEdges(treeEdgesSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type cubeDims(cubeDimsSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type branchLengths(branchLengthsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rootDist(rootDistSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type tipStates(tipStatesSEXP );
        Rcpp::traits::input_parameter< double >::type initial_lambda_01(initial_lambda_01SEXP );
        Rcpp::traits::input_parameter< double >::type initial_lambda_10(initial_lambda_10SEXP );
        Rcpp::traits::input_parameter< double >::type prior_alpha_01(prior_alpha_01SEXP );
        Rcpp::traits::input_parameter< double >::type prior_beta_01(prior_beta_01SEXP );
        Rcpp::traits::input_parameter< double >::type prior_alpha_10(prior_alpha_10SEXP );
        Rcpp::traits::input_parameter< double >::type prior_beta_10(prior_beta_10SEXP );
        Rcpp::traits::input_parameter< int >::type mcmcSize(mcmcSizeSEXP );
        Rcpp::traits::input_parameter< int >::type mcmcBurnin(mcmcBurninSEXP );
        Rcpp::traits::input_parameter< int >::type mcmcSubsample(mcmcSubsampleSEXP );
        NumericVector __result = twoStatePhyloGibbsSampler(treeEdges, cubeDims, branchLengths, rootDist, tipStates, initial_lambda_01, initial_lambda_10, prior_alpha_01, prior_beta_01, prior_alpha_10, prior_beta_10, mcmcSize, mcmcBurnin, mcmcSubsample);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
