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

void modnewton(double *S, double *K, int *nvar, 
	       double *gena, int *nrgen, int *ncgen, double *itmax, double *eps, 
	       double *trace)
{
  int *gset, *cset, ii, jj, kk, ll, nc, ng, ngenElements;
  double *gset2, *gen, *gen2, *gen3, tolin;
  
  ngenElements = *nrgen**ncgen;
  // Rprintf("ngenElements=%d *nvar=%d *itmax=%f\n", ngenElements, *nvar, *itmax);
  
  cset    = (int *)    R_alloc(*nvar,        sizeof(int));
  gset    = (int *)    R_alloc(ngenElements, sizeof(int));
  gset2   = (double *) R_alloc(ngenElements, sizeof(double));
  gen     = (double *) R_alloc(ngenElements, sizeof(double));
  gen2    = (double *) R_alloc(ngenElements, sizeof(double));
  gen3    = (double *) R_alloc(ngenElements, sizeof(double));

  /*   Shift to 0-indexing; */
  for (ii=0;ii<ngenElements;ii++)
    gen[ii] = (gena[ii]-1);
  if(*trace>=2){
    Rprintf("..gen:\n"); printmatd(gen, nrgen, ncgen);
  }
  
  /* Find gset, ng, cset, nc */
  shdunique(gen, &ngenElements, &ng, nvar, gset2);
  for (ii=0; ii<ngenElements;ii++) gset[ii] = (int) gset2[ii]; 
  //Rprintf("gset:"); printveci(gset, &ng);

  complement(gset, &ng, nvar, cset);
  nc = *nvar - ng;

  //Rprintf("cset:"); printveci(cset, &nc);

  double *subt, *S2, *K2, *Sigma2;
  subt   = (double *) R_alloc(ng*ng, sizeof(double));
  S2     = (double *) R_alloc(ng*ng, sizeof(double));
  K2     = (double *) R_alloc(ng*ng, sizeof(double));
  Sigma2 = (double *) R_alloc(ng*ng, sizeof(double));

  double trIS, trISIS;
  double trSS2;
  double Delta2, adj2, prevadj2=0, currparm;
  int itcount=0;
  
  
  // gen2: Version of gen which matches the lower dimensional matrices used later
  for (ii=0;ii<ngenElements;ii++){
    for (jj=0;jj<ng;jj++){
      if (gen[ii]==gset[jj]){
	gen2[ii] = jj;
	gen3[ii] = jj+1;
	break;
      }
    }
  }
  //  Rprintf("gen2:\n"); printmatd(gen2, nrgen, ncgen);
  //  Rprintf("gen3:\n"); printmatd(gen3, nrgen, ncgen);

  // subt <- K(gc)inv(K(cc))K(cc)
  schursubt(K, nvar, gset, &ng, subt);
  // printmatd(subt, &ng, &ng);

  // S2 <- S(gg)
  shdsubmatrix(S, nvar, nvar, gset, &ng, gset, &ng, S2);
  //printmatd(S2, &ng, &ng);

  // trSS2 <- .Call("trAW", gen3, S2)
  trAWprim(gen3, nrgen, ncgen, S2, &ng, &ng, &trSS2);
  //Rprintf("trSS2: %f\n", trSS2);
  
  /*   Iterate here  */

  while(1){
    itcount++;
    // K2 <- K(gg)
    shdsubmatrix(K, nvar, nvar, gset, &ng, gset, &ng, K2);
    //printmatd(K2, &ng, &ng);
    
    // Sigma2    <- solve(K[idx,idx] - subt)
    for (ii=0;ii<ng*ng;ii++){
      Sigma2[ii] = K2[ii]-subt[ii];
    }
    //printmatd(Sigma2, &ng, &ng);
    shdinverse(Sigma2, &ng, &tolin);
    
    trAWprim(gen3, nrgen, ncgen, Sigma2, &ng, &ng, &trIS);
    trAWBWprim(gen3, nrgen, ncgen, Sigma2, &ng, &ng, gen3, nrgen, ncgen, &trISIS);
    //Rprintf("trIS: %f trISIS %f trSS2 %f \n", trIS, trISIS, trSS2);
    
    Delta2  =  trIS - trSS2;
    adj2    =  Delta2 /(trISIS + 0.5*Delta2*Delta2 );

    // Do the update of K
    if (*ncgen==1){
      for (ii=0;ii<*nrgen;ii++){
	kk = (int) gen[ii]*(1+*nvar);
	currparm = K[kk];
	K[kk] = K[kk] + adj2;
      }
    } else {
      for (ii=0;ii<*nrgen;ii++){
	kk = (int) (gen[ii] + *nvar*gen[ii+*nrgen]);
	ll = (int) (gen[ii+*nrgen] + *nvar*gen[ii]);
	currparm = K[kk];
	//Rprintf("gen: (%i %i) kk: %i ll: %i currparm: %f adj: %f\n",
	//(int) gen[ii], (int) gen[ii+*nrgen], kk, ll, currparm, adj2);

	K[kk] = K[kk] + adj2;
	K[ll] = K[ll] + adj2;
      }
    }
    
    if (*trace>=2){
      Rprintf("..Inner Iteration: %i currparm: %15.12f Delta2 %f Update: %15.12f Delta Update: %15.12f\n", 
	      itcount, currparm, Delta2, adj2, abs(adj2-prevadj2));
    }

    //printmatd(K, nvar, nvar);
    if (abs(adj2-prevadj2)<*eps)
      break;
    prevadj2 = adj2;
  }

}


/* Implementation */


SEXP rconipm(SEXP S, SEXP nobs, SEXP K, SEXP Glist, 
	     SEXP maxouter, SEXP maxinner, 
	     SEXP logL, SEXP logLeps, SEXP deltaeps, 
	     SEXP converged,
	     SEXP trace)
{
  R_len_t  nG=length(Glist);
  SEXP ans, Sdims, Kdims,  Gitem, Giidims;
  int  nrS, ncS, nrK, ncK, nrgen, ncgen, ii, jj, Gii, idx;
  double *tracep, *rS, *rK, *rG, *rG2, *ansp, rans;
  double *nobsp, *maxouterp, *maxinnerp, *logLp, *logLepsp, *deltaepsp, *convergedp;
  int outcount=0, outmax;
  double det, trKS, prevlogL=0;
  
  PROTECT(nobs = coerceVector(nobs, REALSXP)) ;
  nobsp = REAL(nobs);

  PROTECT(maxouter = coerceVector(maxouter, REALSXP)) ;
  maxouterp = REAL(maxouter);

  PROTECT(maxinner = coerceVector(maxinner, REALSXP)) ;
  maxinnerp = REAL(maxinner);

  PROTECT(logL = coerceVector(logL, REALSXP)) ;
  logLp = REAL(logL);

  PROTECT(logLeps = coerceVector(logLeps, REALSXP)) ;
  logLepsp = REAL(logLeps);

  PROTECT(deltaeps = coerceVector(deltaeps, REALSXP)) ;
  deltaepsp = REAL(deltaeps);

  PROTECT(converged = coerceVector(converged, REALSXP)) ;
  convergedp = REAL(converged);

  PROTECT(trace = coerceVector(trace, REALSXP)) ;
  tracep = REAL(trace);
  
  Sdims = getAttrib(S, R_DimSymbol);
  if (length(Sdims) < 2) error("Bad Sdims");
  nrS = INTEGER(Sdims)[0];  ncS = INTEGER(Sdims)[1];
  PROTECT(S = AS_NUMERIC(S));
  rS  = REAL(S);

  Kdims = getAttrib(K, R_DimSymbol);
  PROTECT(Kdims = duplicate(Kdims)); 
  if (length(Kdims) < 2) error("Bad Kdims");
  nrK = INTEGER(Kdims)[0];  ncK = INTEGER(Kdims)[1];
  K = duplicate(K);
  PROTECT(K = AS_NUMERIC(K));
  rK  = REAL(K);

  PROTECT(ans = allocVector(REALSXP,nrS*nrS));
  ansp = REAL(ans);

  outmax = (int) *maxouterp;

  while(1){
    outcount++;
    for (Gii=0; Gii<nG; Gii++){
      if (*tracep >= 2)
	Rprintf("..Generator %i\n", Gii);

      PROTECT(Gitem   = AS_NUMERIC(VECTOR_ELT(Glist, Gii)));
      PROTECT(Giidims = getAttrib(Gitem, R_DimSymbol));
      if (length(Giidims) < 2) error("Bad Giidims");
      nrgen   = INTEGER(Giidims)[0]; ncgen   = INTEGER(Giidims)[1];
      rG      = REAL(Gitem);

      modnewton(rS, rK, &nrS, rG, &nrgen, &ncgen,
		maxinnerp, deltaepsp, tracep);
      UNPROTECT(2);
    }

    shddet(rK, &nrS, &det);
    shdtraceAB(rK, &nrS, &nrS, rS, &nrS, &nrS, &trKS);
    *logLp = (*nobsp/2) * (log(det)-trKS); //-(*nobsp/2)* nrS * log(6.283185) + 
    // Rprintf("Outer iteration: %i logL: %f\n", outcount, *logLp);
   
    if (outcount>1){
      if (*tracep >= 1){
	Rprintf("Outer iteration: %i logL: %f diff logL: %f \n",
		outcount, *logLp, *logLp-prevlogL);
      }
      if (*logLp-prevlogL < *logLepsp){
	*convergedp = 1;
	break;
      } else {
	if (outcount==*maxouterp){
	  *convergedp = 0;
	  break;
	}
      }
    }
    prevlogL = *logLp;
  }


  Memcpy(ansp, rK, nrS*nrS);
  setAttrib(ans, R_DimSymbol, Kdims); 
  UNPROTECT(12);  
  return(ans);
}








  //rK[0] = -99999;
  //Memcpy(rK, ansp, nrS*nrS);
  //setAttrib(ans, R_NamesSymbol, getAttrib(K, R_DimSymbol));

