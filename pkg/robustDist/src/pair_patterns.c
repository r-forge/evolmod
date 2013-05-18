#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


SEXP pair_patterns(SEXP x, SEXP y){
  int i, j, nrx,ncx, space_size,*seq_table;
  SEXP ans, Rdimx;
  
  Rdimx = getAttrib(x, R_DimSymbol);
  nrx = INTEGER(Rdimx)[0];
  ncx = INTEGER(Rdimx)[1];

  x = coerceVector(x, INTSXP); 
  y = coerceVector(y, INTSXP); 
  
  seq_table = INTEGER(x);
  space_size = INTEGER(y)[0];
    

  PROTECT(ans = allocMatrix(INTSXP, space_size, space_size));

  for(i=0;i<space_size;i++){
    for(j=0;j<space_size;j++){
      INTEGER(ans)[i + space_size*j]=0; 
    }
  }

  for(i = 0; i < ncx; i++) {
    if ((seq_table[0+nrx*i] !=0) && (seq_table[1+nrx*i] !=0)){
      INTEGER(ans)[seq_table[0+nrx*i]-1 + space_size*(seq_table[1+nrx*i]-1)]++;
    } 
  }
  
  UNPROTECT(1);
  return(ans);
}
