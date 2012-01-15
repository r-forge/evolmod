#include <RcppArmadillo.h>

using namespace Rcpp; 

 /* Arguments: 
                tE: edge matrix of the ape object phylo
                nIN: number of internal nodes (not necessary, but convinient to get it from phylo)
                tS: integer vector of tip states (-1=missing value)
                nS: state space size (e.g. 2 for a binary trait)
                pM: array of probability matrices for each edge of the tree

     Two important assumptions:
                1. edges in the the edge matrix and probability matrices are in the "pruningwise" order; 
                   see ?reorder.phylo for more details
                2. tip state vector is ordered according to the tip numbering in the edge matrix
  */

arma::mat part_likelihoods(IntegerMatrix treeEdges, 
			int numIntNodes, 
			IntegerVector tipStates,
			int numStates,
			arma::cube cubeProbMat
			){

  //IntegerMatrix treeEdges(tE);

  // get the number of edges
  int numEdges = treeEdges.nrow();
  //int numIntNodes = as<int>(nIN);
  //IntegerVector tipStates(tS);
  //int numStates = as<int>(nS);
  //NumericVector vecProbMat(pM);
  //arma::cube cubeProbMat(vecProbMat.begin(), numStates, numStates, numEdges, false);
    
  // get the number of tips in the tree
  int numTips = tipStates.size();
  
  // prepare a matrix for storing regular (backward) partial likelihoods
  arma::mat partialLike = arma::zeros<arma::mat>(numTips + numIntNodes, numStates);

  for (int i=0; i < numTips; i++){
    if (tipStates[i] == -1){// -1 denotes a missing value
      partialLike.row(i) = arma::ones<arma::rowvec>(numStates);
    }else{
      partialLike(i, tipStates[i]) = 1.0;
    }
  }

  // compute regular partial likelihoods for all internal nodes
  for (int i=0; i < numEdges; i+=2){      
    // parent1 = treeEdges[i,0] or treeEdges[i+1,0] also treeEdges indices should be shifted down by one
    partialLike.row(treeEdges(i,0)-1) = (partialLike.row(treeEdges(i,1)-1)*cubeProbMat.slice(i).t())%(partialLike.row(treeEdges(i+1,1)-1)*cubeProbMat.slice(i+1).t());            
  }

  return partialLike;
}
