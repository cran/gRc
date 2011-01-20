#include <Rdefines.h>
#include <R.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "shd_print.h"
#include "shd_utils.h"
#include "clevertrace.h"

/* Forward declarations */

#define	abs(x)			((x) >= 0 ? (x) : -(x))

SEXP scorematching(
		   // Empirical covariance matrix
		   SEXP S, 
		   // Number of cases
		   SEXP nobs, 
		   // Generator list
		   SEXP vcc,
		   SEXP ecc,
		   // log likelihood
		   SEXP logL, 
		   // debugginf info
		   SEXP debug
		   )
{
  R_len_t  nvcc=length(vcc);
  R_len_t  necc=length(ecc);
  int nG = nvcc + necc;
  
  SEXP result, Sdims;
  SEXP Uitem, Uiidims, Vitem, Viidims;
  int nrU, ncU, nrV, ncV;
  double *rU, *rV;

  int     nrS, ncS,  ii;
  double  *rS,  *resultp;
  double  *nobsp, *logLp,  *debugp;
  
  double *Amat, *Bvec, tmp;
  int uu,vv; 

  PROTECT(Sdims = getAttrib(S, R_DimSymbol));
  if (length(Sdims) < 2) error("Bad Sdims");
  nrS = INTEGER(Sdims)[0];  
  ncS = INTEGER(Sdims)[1];
  
  PROTECT(S = AS_NUMERIC(S));  
  //S   = duplicate(S);
  rS  = REAL(S);

  PROTECT(nobs = coerceVector(nobs, REALSXP)) ;
  nobsp = REAL(nobs);

  PROTECT(logL = coerceVector(logL, REALSXP)) ;
  logLp = REAL(logL);

  PROTECT(debug = coerceVector(debug, REALSXP)) ;
  debugp = REAL(debug);

  PROTECT(result = allocVector(REALSXP,nG));
  resultp = REAL(result);
  /*   Memcpy(resultp, rK, nrS*nrS); */

  //outmax = (int) *maxouterp;

  Amat  = (double *) R_alloc(nG*nG, sizeof(double)); 
  Bvec  = (double *) R_alloc(nG,    sizeof(double)); 
  for (ii=0;ii<nG*nG;ii++)
    Amat[ii]=0;
  for (ii=0;ii<nG;ii++)
    Bvec[ii]=0;
  //printmatd(Amat, &nG, &nG);

  // VCC x VCC
  for (uu=0; uu<nvcc; uu++){
    PROTECT(Uitem   = AS_NUMERIC(VECTOR_ELT(vcc, uu)));
    PROTECT(Uiidims = getAttrib(Uitem, R_DimSymbol));
    if (length(Uiidims) < 2) error("Bad Uiidims");
    nrU   = INTEGER(Uiidims)[0]; 
    ncU   = INTEGER(Uiidims)[1];
    rU    = REAL(Uitem);

    trAWBprim(rU, &nrU, &ncU, rS, &nrS, &ncS, rU, &nrU, &ncU, &tmp); 
    Amat[uu*(nG+1)] = tmp; // Upper left
    Bvec[uu]        = nrU;
    //printmatd(rG, &nrgen, &ncgen);
    //Rprintf("value: %f\n", tmp); 
    UNPROTECT(2);
  }

  //printmatd(Amat, &nG, &nG);

  if (necc>0){
    // VCC x ECC
    for (uu=0; uu<nvcc; uu++){
      PROTECT(Uitem   = AS_NUMERIC(VECTOR_ELT(vcc, uu)));
      PROTECT(Uiidims = getAttrib(Uitem, R_DimSymbol));
      if (length(Uiidims) < 2) error("Bad Uiidims");
      nrU   = INTEGER(Uiidims)[0]; 
      ncU   = INTEGER(Uiidims)[1];
      rU    = REAL(Uitem);
      // Rprintf("Uterm:\n"); printmatd(rU, &nrU, &ncU);
      for (vv=0; vv<necc; vv++){
	
	PROTECT(Vitem   = AS_NUMERIC(VECTOR_ELT(ecc, vv)));
	PROTECT(Viidims = getAttrib(Vitem, R_DimSymbol));
	if (length(Viidims) < 2) error("Bad Viidims");
	nrV   = INTEGER(Viidims)[0]; 
	ncV   = INTEGER(Viidims)[1];
	rV    = REAL(Vitem);   
	// Rprintf("Vterm:\n"); printmatd(rV, &nrV, &ncV);
	
	trAWBprim(rU, &nrU, &ncU, rS, &nrS, &ncS, rV, &nrV, &ncV, &tmp); 
	Amat[uu + nG*(vv+nvcc)] = tmp; // Upper right      
	Amat[vv+nvcc + nG*uu]   = tmp; // Lower left
	//Rprintf("value: %f\n", tmp); 
	UNPROTECT(2);
      }
      UNPROTECT(2);
    }
    //printmatd(Amat, &nG, &nG);
    
    
    // ECC x ECC
    for (uu=0; uu<necc; uu++){
      PROTECT(Uitem   = AS_NUMERIC(VECTOR_ELT(ecc, uu)));
      PROTECT(Uiidims = getAttrib(Uitem, R_DimSymbol));
      if (length(Uiidims) < 2) error("Bad Uiidims");
      nrU   = INTEGER(Uiidims)[0]; 
      ncU   = INTEGER(Uiidims)[1];
      rU    = REAL(Uitem);
      //Rprintf("Uterm:\n"); printmatd(rU, &nrU, &ncU);
      for (vv=uu; vv<necc; vv++){
	
	PROTECT(Vitem   = AS_NUMERIC(VECTOR_ELT(ecc, vv)));
	PROTECT(Viidims = getAttrib(Vitem, R_DimSymbol));
	if (length(Viidims) < 2) error("Bad Viidims");
	nrV   = INTEGER(Viidims)[0]; 
	ncV   = INTEGER(Viidims)[1];
	rV    = REAL(Vitem);   
	//Rprintf("Vterm:\n"); printmatd(rV, &nrV, &ncV);
	
	trAWBprim(rU, &nrU, &ncU, rS, &nrS, &ncS, rV, &nrV, &ncV, &tmp); 
	Amat[uu+nvcc + nG*(vv+nvcc)] = tmp; // Upper right of lower right      
	Amat[vv+nvcc + nG*(uu+nvcc)]   = tmp; // Lower left of lower right
	//Rprintf("value: %f\n", tmp); 
	UNPROTECT(2);
      }
      UNPROTECT(2);
    }
  }
  //printmatd(Amat, &nG, &nG);
  //printvecd(Bvec, &nG);

  int ncB=1;
  double tolin=0.00001;
  shdsolve(Amat, &nG, Bvec, &ncB, &tolin);

  //Rprintf("ANSWER:\n");printvecd(Bvec, &nG);


  //  setAttrib(result, R_DimSymbol, Kdims); 
  UNPROTECT(6);
  //Rprintf("Exiting rconipm...\n");  
  Memcpy(resultp, Bvec, nG);
  return(result);
}



