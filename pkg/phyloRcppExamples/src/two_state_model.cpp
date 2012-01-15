#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat inf_generator(double lambda_01, double lambda_10){

  arma::mat generator = arma::mat(2,2);
  generator(0,0) = -lambda_01;
  generator(0,1) = lambda_01;
  generator(1,0) = lambda_10;
  generator(1,1) = -lambda_10;
  
  return generator;
}

//Computes finite time transition probabilities using analytic formulae.

arma::mat trans_prob(double lambda_01, double lambda_10, double time){

  double total_rate = lambda_01 + lambda_10;               

  arma::mat prob_matrix = arma::mat(2,2);
  
  prob_matrix(0,0) = (lambda_10 + lambda_01*exp(-total_rate*time))/total_rate;
  prob_matrix(0,1) = (lambda_01 - lambda_01*exp(-total_rate*time))/total_rate;
  prob_matrix(1,0) = (lambda_10 - lambda_10*exp(-total_rate*time))/total_rate;
  prob_matrix(1,1) = (lambda_01 + lambda_10*exp(-total_rate*time))/total_rate;

  return prob_matrix;
}

/* 
  Computes and eigen decomposition of the CTMC
  infinitesimal generator and saves the
  eigenvalues, left eigenvectors (in matrix form)
  and the inverse of the above matrix. 
  Uses analytic formulae for the two state model.
*/

List eigen_decomp(double lambda_01, double lambda_10){
                   
  double total_rate = lambda_01 + lambda_10;
  NumericVector eigen_values = NumericVector::create(0, -total_rate);
                   
  arma::mat left_eigen = arma::mat(2,2);
  arma::mat right_eigen = arma::mat(2,2);
                   
  right_eigen(0,0) = 1;
  right_eigen(1,0) = 1;
  right_eigen(0,1) = lambda_01/total_rate;
  right_eigen(1,1) = -lambda_10/total_rate;

  left_eigen(0,0) = lambda_10/total_rate;
  left_eigen(0,1) = lambda_01/total_rate;
  left_eigen(1,0) = 1;
  left_eigen(1,1) = -1;

  List eigen_decomp;
  
  eigen_decomp["evalues"] = eigen_values;
  eigen_decomp["evectors"] = right_eigen;
  eigen_decomp["einvvectors"] = left_eigen;
  
  return eigen_decomp;
}


// The functions below expose the above internal C++ functions to R.
// They are here only for testing purposes.

RcppExport SEXP inf_generator_to_R(SEXP lambda_01, SEXP lambda_10){
  return wrap(inf_generator(as<double>(lambda_01), as<double>(lambda_10)));
}

RcppExport SEXP trans_prob_to_R(SEXP lambda_01, SEXP lambda_10, SEXP time){
  return wrap(trans_prob(as<double>(lambda_01), as<double>(lambda_10), as<double>(time)));
}

RcppExport SEXP eigen_decomp_to_R(SEXP lambda_01, SEXP lambda_10, SEXP time){
  return wrap(eigen_decomp(as<double>(lambda_01), as<double>(lambda_10)));
}

RCPP_MODULE(two_state){ 
  function( "inf_generator", &inf_generator );
  function( "trans_prob", &trans_prob );
  function( "eigen_decomp", &eigen_decomp );
}
