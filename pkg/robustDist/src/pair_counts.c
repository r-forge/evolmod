#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


SEXP pair_counts(SEXP obj1, SEXP obj2, SEXP obj3, SEXP obj4){
  int i, j, nrx, seq_length, space_size, *seq_table;
  double *cond_counts, *marg_counts, stat_counts, total_counts;
  SEXP ans, Rdimseq, Rdimspace;
  
  Rdimseq = getAttrib(obj1, R_DimSymbol);
  nrx = INTEGER(Rdimseq)[0];
  seq_length = INTEGER(Rdimseq)[1];

  if (nrx != 2) error("The number of sequences must be 2!");

  Rdimspace = getAttrib(obj2, R_DimSymbol);
  space_size = INTEGER(Rdimspace)[1];

  obj1 = coerceVector(obj1, INTSXP); 
  obj2 = coerceVector(obj2, REALSXP); 
  obj3 = coerceVector(obj3, REALSXP); 
  obj4 = coerceVector(obj4, REALSXP); 
  
  seq_table = INTEGER(obj1);
  cond_counts = REAL(obj2);
  marg_counts = REAL(obj3);
  stat_counts = REAL(obj4)[0];

    

  PROTECT(ans = allocVector(REALSXP,1));

  REAL(ans)[0]=0;

  for(i = 0; i < seq_length; i++) {
    if ((seq_table[0+2*i] !=0) && (seq_table[1+2*i] !=0)){
      REAL(ans)[0] += cond_counts[seq_table[0+2*i]-1 + space_size*(seq_table[1+2*i]-1)];
    }else{
      if ((seq_table[0+2*i] ==0) && (seq_table[1+2*i] !=0)){
	REAL(ans)[0] += marg_counts[seq_table[1+2*i]-1];
      }else{
	if ((seq_table[0+2*i] !=0) && (seq_table[1+2*i] ==0)){
	  REAL(ans)[0] += marg_counts[seq_table[0+2*i]-1];
	}else{
	  REAL(ans)[0] += stat_counts;
	}
      } 
    }
  }
  
  UNPROTECT(1);
  return(ans);
}
