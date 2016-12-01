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

void modnewton(double *S, 
	       double *K, 
	       int *nvar, 
	       double *gena, int *nrgen, int *ncgen, 
	       // Maximum iterations allowed
	       double *itmax, 
	       // Iterations performed
	       int *inner, 
	       // Convergence criterion; |parmchange|
	       double *eps, 
	       // |parmchange|
	       double *parmchg,
	       double *trace,
	       int *gcsetI, 
	       double *gcsetD,
	       double *NRtmp,
	       double *WORK,
	       double *GSETtmp	       
	       )
{
  int *gset, *cset, ii, jj, kk, ll, nc, ng, ngenElements;
  double *gset2, *gen, *gen2, *gen3, tolin;
  double *subt, *S2, *K2, *Sigma2;


  ngenElements = *nrgen**ncgen; 
  // Rprintf("ngenElements=%d *nvar=%d *itmax=%f\n", ngenElements, *nvar, *itmax);

  gen     = GSETtmp;
  gen2    = GSETtmp + 1*ngenElements;
  gen3    = GSETtmp + 2*ngenElements;
  gset2   = GSETtmp + 3*ngenElements;


  /*   Shift to 0-indexing; */
  for (ii=0;ii<ngenElements;ii++)
    gen[ii] = (gena[ii]-1);
  if(*trace>=6){
    Rprintf("......gen:\n"); printmatd(gen, nrgen, ncgen);
  }

  /* Find gset, ng, cset, nc */


  shdunique(gen, &ngenElements, &ng, nvar, gcsetD);
  //goto EXIT;        
  //Rprintf("unique ok...\n");
  for (ii=0; ii<*nvar;ii++) 
    gcsetI[ii] = (int) gcsetD[ii]; 

  gset = gcsetI;
  cset = gcsetI + ng;
  nc   = *nvar - ng;
  
  //Rprintf("gset:\n"); printveci(gset, &ng);
  //Rprintf("cset:\n"); printveci(cset, &nc);

  subt   = NRtmp;
  S2     = NRtmp + 1*ng*ng;
  K2     = NRtmp + 2*ng*ng;
  Sigma2 = NRtmp + 3*ng*ng;
  
  double trIS, trISIS, trSS2;
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
  //Rprintf("Done creating gen2,gen3\n");

  // subt <- K(gc)inv(K(cc))K(cc)
  schursubt2(K, nvar, gset, &ng, cset, &nc, WORK, subt);
  //Rprintf("schursubt2 ok...\n");
  // printmatd(subt, &ng, &ng);

  // S2 <- S(gg)
  shdsubmatrix(S, nvar, nvar, gset, &ng, gset, &ng, S2);
  //Rprintf("S2 ok...\n");
  //printmatd(S2, &ng, &ng);

  // trSS2 <- .Call("trAW", gen3, S2)
  trAWprim(gen3, nrgen, ncgen, S2, &ng, &ng, &trSS2);
  //Rprintf("trAW ok...\n");
  //Rprintf("trSS2: %f\n", trSS2);
  
  /* ***  Iterate here *** */

  while(1){ 
    //Rprintf("while...\n");
    itcount++;
    // K2 <- K(gg)
    shdsubmatrix(K, nvar, nvar, gset, &ng, gset, &ng, K2);
    //Rprintf("K2 ok...\n");printmatd(K2, &ng, &ng);
    
    // Sigma2    <- solve(K[idx,idx] - subt)
    for (ii=0;ii<ng*ng;ii++){
      Sigma2[ii] = K2[ii]-subt[ii];
    }
    //Rprintf("Sigma2 ok...\n"); printmatd(Sigma2, &ng, &ng);

    shdinverse(Sigma2, &ng, &tolin);
    //Rprintf("inv(Sigma2) ok...\n"); printmatd(Sigma2, &ng, &ng);

    //printmatd(gen3, nrgen, ncgen);
    trAWprim(gen3, nrgen, ncgen, Sigma2, &ng, &ng, &trIS);
    trAWBWprim(gen3, nrgen, ncgen, Sigma2, &ng, &ng, gen3, nrgen, ncgen, &trISIS);
    //if (*trace>=4)
    //Rprintf("....trIS: %f trISIS %f trSS2 %f \n", trIS, trISIS, trSS2);
    
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
    
    if (*trace>=6){
      Rprintf("......Inner Iteration: %i currparm: %15.12f Delta2 %f Update: %15.12f Delta Update: %15.12f\n",
	      itcount, currparm, Delta2, adj2, abs(adj2-prevadj2));
    }

    *parmchg = abs(adj2-prevadj2);
    
    if (*parmchg < *eps){ 
      break;
    }
    prevadj2 = adj2;
  } // while 
  *inner = itcount;
  
/*   if (*trace>=3){ */
/*     Rprintf("...Inner loop done; iterations: %i \n", itcount); */
/*   } */

  // EXIT: Rprintf("OUT...\n");
}


/* Implementation */


SEXP rconipm(SEXP S, SEXP nobs, SEXP K, SEXP Glist, 
	     SEXP maxouter, SEXP maxinner, 
	     SEXP logLeps, SEXP deltaeps, 
	     SEXP converged,
	     SEXP debug)
{
  R_len_t  nG=length(Glist);
  SEXP Kans, Sdims, Kdims,  Gitem, Giidims;
  SEXP outer, globalinner;
  SEXP logL;
  double *outerp, *globalinnerp;
  int     nrS, ncS, nrK, ncK, nrgen, ncgen,  Gii;
  double  *rS, *rK, *rG, *Kansp;
  double  *nobsp, *maxouterp, *maxinnerp, *logLp, *logLepsp, *deltaepsp, *convergedp, *debugp, parmchg;
  int     outcount=0, outmax, inner;
  double  det, trKS, prevlogL=0;
  
  int ngenElements, ng, nc, maxng=0, maxnc=0, maxngenElements=0;
  double *gcsetD, *NRtmp, *WORK, *GSETtmp;
  int    *gcsetI;
  

  PROTECT(nobs = coerceVector(nobs, REALSXP)) ;
  nobsp = REAL(nobs);

  PROTECT(maxouter = coerceVector(maxouter, REALSXP)) ;
  maxouterp = REAL(maxouter);


  PROTECT(maxinner = coerceVector(maxinner, REALSXP)) ;
  maxinnerp = REAL(maxinner);

  PROTECT(logL = allocVector(REALSXP, 1)) ;
  logLp = REAL(logL);

  PROTECT(logLeps = coerceVector(logLeps, REALSXP)) ;
  logLepsp = REAL(logLeps);

  PROTECT(deltaeps = coerceVector(deltaeps, REALSXP)) ;
  deltaepsp = REAL(deltaeps);

  PROTECT(converged = coerceVector(converged, REALSXP)) ;
  convergedp = REAL(converged);

  PROTECT(debug = coerceVector(debug, REALSXP)) ;
  debugp = REAL(debug);

  /* *** S *** */
  PROTECT(Sdims = getAttrib(S, R_DimSymbol));
  if (length(Sdims) < 2) error("Bad Sdims");
  nrS = INTEGER(Sdims)[0];  
  ncS = INTEGER(Sdims)[1];

  //S   = duplicate(S);  
  PROTECT(S = AS_NUMERIC(S));  
  rS  = REAL(S);
  
  /* *** K *** */
  PROTECT(Kdims = getAttrib(K, R_DimSymbol));
  if (length(Kdims) < 2) error("Bad Kdims");
  nrK = INTEGER(Kdims)[0];  
  ncK = INTEGER(Kdims)[1];

  //K   = duplicate(K); 
  PROTECT(K = AS_NUMERIC(K));
  rK  = REAL(K);

  PROTECT(Kans = allocVector(REALSXP,nrS*nrS));
  Kansp = REAL(Kans);
  Memcpy(Kansp, rK, nrS*nrS);

  PROTECT(outer=allocVector(REALSXP, 1));
  outerp = REAL(outer);

  PROTECT(globalinner=allocVector(REALSXP, 1));
  globalinnerp = REAL(globalinner);

  outmax = (int) *maxouterp;

  // Allocate the necessary space for temporary variables 
  // 
  int nMAX = (int) nrK*(nrK-1)/2;
  gcsetI  = (int *)    R_alloc(nMAX, sizeof(int));
  gcsetD  = (double *) R_alloc(nMAX, sizeof(double));

  //gset    = (double *) malloc(nrK*sizeof(double));
  for (Gii=0; Gii<nG; Gii++){
    
    PROTECT(Gitem   = AS_NUMERIC(VECTOR_ELT(Glist, Gii)));
    PROTECT(Giidims = getAttrib(Gitem, R_DimSymbol));
    if (length(Giidims) < 2) error("Bad Giidims");
    nrgen   = INTEGER(Giidims)[0]; 
    ncgen   = INTEGER(Giidims)[1];    
    rG      = REAL(Gitem);
    
    ngenElements = nrgen * ncgen;
    shdunique(rG, &ngenElements, &ng, &nMAX, gcsetD);
    nc = nrK - ng;
    
    if (ng>maxng)
      maxng = ng;
    if (nc>maxnc)
      maxnc = nc;
    if (ngenElements>maxngenElements)
      maxngenElements = ngenElements;
    
    UNPROTECT(2);
  }

  //Rprintf("nrK: %i maxngenElements: %i, maxng; %i maxnc: %i\n",
  //  nrK, maxngenElements, maxng, maxnc);
  
  NRtmp   = (double *) R_alloc(4*maxng*maxng, sizeof(double));
  WORK    = (double *) R_alloc(nrK*nrK, sizeof(double));
  GSETtmp = (double *) R_alloc(4*maxngenElements, sizeof(double));    
    
  //
  // Memory allocation done... !!
  
  if (*debugp>=1){
    shddet(rK, &nrS, &det);    //Rprintf("det %f\n", det);
    shdtraceAB(rK, &nrS, &nrS, rS, &nrS, &nrS, &trKS); // Rprintf("trKS %f\n", trKS);
    *logLp = (*nobsp/2) * (log(det)-trKS); //-(*nobsp/2)* nrS * log(6.283185) + 
    Rprintf(".Before iteration: %i logL: %f\n", outcount, *logLp);
    prevlogL = *logLp;
  }

  while(1){
    outcount++;
    globalinner = 0; // The maximal # of inner iterations.
    for (Gii=0; Gii<nG; Gii++){
/*       if (*debugp >= 2) */
/* 	Rprintf("..Generator %i\n", Gii); */

      PROTECT(Gitem   = AS_NUMERIC(VECTOR_ELT(Glist, Gii)));
      PROTECT(Giidims = getAttrib(Gitem, R_DimSymbol));
      if (length(Giidims) < 2) error("Bad Giidims");
      nrgen   = INTEGER(Giidims)[0]; 
      ncgen   = INTEGER(Giidims)[1];
      rG      = REAL(Gitem);

      modnewton(rS, Kansp, &nrS, rG, &nrgen, &ncgen, maxinnerp, &inner, 
 		deltaepsp, 
		&parmchg,
		debugp, 
 		gcsetI, gcsetD, NRtmp, WORK, GSETtmp); 
      UNPROTECT(2);

      //printmatd(rG, &nrgen, &ncgen);
      if (*debugp>=4){
	shddet(Kansp, &nrS, &det);    //Rprintf("det %f\n", det);
	shdtraceAB(Kansp, &nrS, &nrS, rS, &nrS, &nrS, &trKS); // Rprintf("trKS %f\n", trKS);
	*logLp = (*nobsp/2) * (log(det)-trKS); //-(*nobsp/2)* nrS * log(6.283185) + 
	Rprintf("....Generator %3i Inner iterations %i logL: %f\n", Gii, inner, *logLp);
      }
      if (inner>*globalinnerp) // The maximal # of inner iterations.
	*globalinnerp = inner;
    }

    shddet(Kansp, &nrS, &det);    //Rprintf("det %f\n", det);
    shdtraceAB(Kansp, &nrS, &nrS, rS, &nrS, &nrS, &trKS); // Rprintf("trKS %f\n", trKS);
    *logLp = (*nobsp/2) * (log(det)-trKS); //-(*nobsp/2)* nrS * log(6.283185) + 
   
    if (*debugp >= 1){
      Rprintf(".Outer iteration: %3i logL: %20.12f diff logL: %18.12f \n",
		outcount, *logLp, *logLp-prevlogL);
    } 

    if (outcount>1){
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

  *outerp = outcount;
  //printmatd(Kansp, &nrS, &nrS);
  setAttrib(Kans, R_DimSymbol, Kdims); 

  // NB: only elements 0 and 3 are used currently
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, Kans);           // K
  SET_VECTOR_ELT(res, 1, outer);          // # of outer iterations
  SET_VECTOR_ELT(res, 2, globalinner);    // max # of inner iterations
  SET_VECTOR_ELT(res, 3, logL);
  UNPROTECT(16);
  return(res);
}




  // Free(gcsetD); Free(gcsetI); Free(NRtmp); Free(WORK); Free(GSETtmp);



/*   PROTECT(S = AS_NUMERIC(S)); */
/*   Sdims = getAttrib(S, R_DimSymbol); */
/*   PROTECT(Sdims = duplicate(Sdims));  */
/*   if (length(Sdims) < 2) error("Bad Sdims"); */
/*   nrS = INTEGER(Sdims)[0];  ncS = INTEGER(Sdims)[1];   */
/*   //  S   = duplicate(S); */

/*   rS  = REAL(S); */
  
/*   /\* *** K *** *\/ */
/*   PROTECT(K = AS_NUMERIC(K)); */
/*   Kdims = getAttrib(K, R_DimSymbol); */
/*   PROTECT(Kdims = duplicate(Kdims));  */
/*   if (length(Kdims) < 2) error("Bad Kdims"); */
/*   nrK = INTEGER(Kdims)[0];  ncK = INTEGER(Kdims)[1];   */
/*   //K   = duplicate(K); */
/*   rK  = REAL(K); */


