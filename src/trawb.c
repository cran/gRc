#include <Rdefines.h>
#include <R.h>

SEXP trawb(SEXP A, SEXP W, SEXP B)
{
  int nrA, ncA, nrW, ncW, nrB, ncB;
  /*   int k1, k2; */
  int i, j;
  double *rans, *rA=REAL(A), *rW=REAL(W), *rB=REAL(B);
  SEXP Adims, Wdims, Bdims, ans;

  /* int a,b,g,d; */
  int aa,bb,gg,dd;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);
  Bdims = getAttrib(B, R_DimSymbol);

  nrA = INTEGER(Adims)[0];
  ncA = INTEGER(Adims)[1];
  
  nrW = INTEGER(Wdims)[0];
  ncW = INTEGER(Wdims)[1];

  nrB = INTEGER(Bdims)[0];
  ncB = INTEGER(Bdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 0.0;

  if (ncA==2){
    if (ncB==2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	bb = (int)  rA[i + nrA*1]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  dd = (int)  rB[j + nrB*1]-1;
	  
	  *rans = *rans +
	    rW[(int) (bb+nrW*gg)]*(aa==dd) +
	    rW[(int) (aa+nrW*gg)]*(bb==dd) +
	    rW[(int) (bb+nrW*dd)]*(aa==gg) +
	    rW[(int) (aa+nrW*dd)]*(bb==gg) ;
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	bb = (int)  rA[i + nrA*1]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;

	  *rans = *rans +
	    rW[(int) (gg+nrW*aa)]*(gg==bb) +
	    rW[(int) (gg+nrW*bb)]*(gg==aa);
	}
      ;
      }
    }
  } else { /* ncA==1 */
    if (ncB==2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  dd = (int)  rB[j + nrB*1]-1;
	  
	  *rans = *rans +
	    rW[(int) (aa+nrW*gg)]*(aa==dd) +
	    rW[(int) (aa+nrW*dd)]*(aa==gg);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  
	  *rans = *rans +
	    rW[(int) (aa+nrW*aa)]*(aa==gg);
	}
      }

    }
  }
  
  UNPROTECT(1);
  return(ans);
}


SEXP trawbw(SEXP A, SEXP W, SEXP B)
{
  int nrA, ncA, nrW, ncW, nrB, ncB;
  /*   int k1, k2; */
  int i, j;
  double *rans, *rA=REAL(A), *rW=REAL(W), *rB=REAL(B);
  SEXP Adims, Wdims, Bdims, ans;

  /*   int a,b,g,d; */
  int aa,bb,gg,dd;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);
  Bdims = getAttrib(B, R_DimSymbol);

  nrA = INTEGER(Adims)[0];
  ncA = INTEGER(Adims)[1];
  
  nrW = INTEGER(Wdims)[0];
  ncW = INTEGER(Wdims)[1];

  nrB = INTEGER(Bdims)[0];
  ncB = INTEGER(Bdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 0.0;

  if (ncA==2){
    if (ncB==2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	bb = (int)  rA[i + nrA*1]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  dd = (int)  rB[j + nrB*1]-1;
	  
	  *rans = *rans +	
	    2*(rW[(int) (bb+nrW*gg)]*rW[(int) (aa+nrW*dd)] +
	    rW[(int) (aa+nrW*gg)]*rW[(int) (bb+nrW*dd)]);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	bb = (int)  rA[i + nrA*1]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;

	  *rans = *rans +	
	    2*(rW[(int) (aa+nrW*gg)]*rW[(int) (bb+nrW*gg)]);
	}
      ;
      }
    }
  } else { /* ncA==1 */
    if (ncB==2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  dd = (int)  rB[j + nrB*1]-1;
	  
	  *rans = *rans +	
	    2*(rW[(int) (aa+nrW*gg)]* rW[(int) (aa+nrW*dd)]);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i + nrA*0]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j + nrB*0]-1;
	  
	  *rans = *rans +	
	    rW[(int) (aa+nrW*gg)]*rW[(int) (aa+nrW*gg)];
	}
      }
    }
  }
  
  UNPROTECT(1);
  return(ans);
}




  /*   Rprintf(" nrA %d ncA %d\n nrW %d ncW %d\n nrB %d ncB %d \n", */
  /* 	  nrA, ncA, nrW, ncW, nrB, ncB); */

  /*   int Amat[nrA][ncA]; */
  /*   for (i=0; i<nrA; i++){ */
  /*     for (j=0; j<ncA; j++){ */
  /*       Amat[i][j] = rA[i + nrA*j]; */
  /*       /\*       Rprintf("A: %i-%f", Amat[i][j], rA[i + nrA*j]);  *\/ */
  /*     } */
  /*     /\*     Rprintf("\n");  *\/ */
  /*   } */
  
  /*   /\*   Rprint("B:\n"); *\/ */
  /*   int Bmat[nrB][ncB]; */
  /*   for (i=0; i<nrB; i++){ */
  /*     for (j=0; j<ncB; j++){ */
  /*       Bmat[i][j] = rB[i + nrB*j]; */
  /*       /\*       Rprintf(" %i", Bmat[i][j]); *\/ */
  /*     } */
  /*     /\*     Rprintf("\n"); *\/ */
  /*   } */


	/*       g  = (int) (Bmat[j][0]-1); */
	/*       d  = (int) (Bmat[j][1]-1); */
	/*       Rprintf("B entries: %i %i \n", Bmat[j][0], Bmat[j][1]);  */
	/*       Rprintf("B entries: %i %i %i %i\n", g,d, gg, dd); */
	/*       Rprintf("%i %i \n", (b+nrW*g),  (bb+nrW*gg) ); */


      /*     a  = (int) (Amat[i][0]-1); */
      /*     b  = (int) (Amat[i][1]-1); */
      /*     Rprintf("A entries: %i %i\n", Amat[i][0], Amat[i][1]); */
      /*     Rprintf("A entries: %i %i %i %i\n", a,b, aa, bb); */
      /*     k1 = (Amat[i][0]-1) + nrW*(Amat[i][1]-1); */




/* 	rW[((int) (bb+nrW*gg))]*(a==d) + */
/* 	rW[((int) (aa+nrW*gg))]*(b==d) + */
/* 	rW[((int) (bb+nrW*dd))]*(a==g) + */
/* 	rW[((int) (aa+nrW*dd))]*(b==g) ; */
/* 	rW[((int) (b+nrW*g))]*(a==d) + */
/* 	rW[((int) (a+nrW*g))]*(b==d) + */
/* 	rW[((int) (b+nrW*d))]*(a==g) + */
/* 	rW[((int) (a+nrW*d))]*(b==g) ; */

/*       Rprintf("kk %d w %f\n", kk, xxx); */
/*       k2 = (Bmat[j][0]-1) + nrW*(Bmat[j][1]-1); */
/*       Rprintf("k1 %d k2 %d\n", k1, k2); */
