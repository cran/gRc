#include <Rdefines.h>
#include <R.h>

void trAWprim(double *rA, int *nrA, int *ncA,
	      double *rW, int *nrW, int *ncW, double *rans);

void trAWBprim(double *rA, int *nrA, int *ncA,
	       double *rW, int *nrW, int *ncW,
	       double *rB, int *nrB, int *ncB, double *rans);

void trAWBWprim(double *rA, int *nrA, int *ncA,
		double *rW, int *nrW, int *ncW,
		double *rB, int *nrB, int *ncB, double *rans);

void trAWBVprim(double *rA, int *nrA, int *ncA,
		double *rW, int *nrW, int *ncW,
		double *rB, int *nrB, int *ncB, 
		double *rV, int *nrV, int *ncV, double *rans);


/* Implementation */

void trAWprim(double *rA, int *nrA, int *ncA,
	      double *rW, int *nrW, int *ncW, double *rans)
{
  int ii,idx;
  *rans = 0;
  if (*ncA==2){
    for (ii=0; ii<*nrA; ii++){
      idx = (int) (rA[ii]-1 + *nrW* (rA[ii+*nrA]-1));
      //Rprintf("%i ", idx);
      *rans = *rans + rW[idx];
    }
    *rans = 2* *rans;    
  } else {
    for (ii=0; ii<*nrA; ii++){
      idx = (int) ((rA[ii]-1)*(*nrW+1));
      //Rprintf("%i %f \n", idx, rA[ii]);
      *rans = *rans + rW[idx];
    }
  }
}


void trAWBprim(double *rA, int *nrA, int *ncA,
	       double *rW, int *nrW, int *ncW,
	       double *rB, int *nrB, int *ncB, double *rans)
{
  int i,j;
  int aa,bb,gg,dd;
  *rans=0;
  
  //printmatd(rA, nrA, ncA);
  //printmatd(rB, nrB, ncB);
  
  if (*ncA==2){
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;
	  
	  *rans = *rans +
	    rW[(int) (bb+*nrW*gg)]*(aa==dd) +
	    rW[(int) (aa+*nrW*gg)]*(bb==dd) +
	    rW[(int) (bb+*nrW*dd)]*(aa==gg) +
	    rW[(int) (aa+*nrW*dd)]*(bb==gg) ;
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;

	  *rans = *rans +
	    rW[(int) (gg+*nrW*aa)]*(gg==bb) +
	    rW[(int) (gg+*nrW*bb)]*(gg==aa);
	}
      ;
      }
    }
  } else { /* ncA==1 */
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;
	  
	  *rans = *rans +
	    rW[(int) (aa+*nrW*gg)]*(aa==dd) +
	    rW[(int) (aa+*nrW*dd)]*(aa==gg);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  
	  *rans = *rans +
	    rW[(int) (aa+*nrW*aa)]*(aa==gg);
	}
      }

    }
  }

}


void trAWBWprim(double *rA, int *nrA, int *ncA,
	       double *rW, int *nrW, int *ncW,
	       double *rB, int *nrB, int *ncB, double *rans)
{

  int i,j;
  int aa,bb,gg,dd;
  *rans=0;

  if (*ncA==2){
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;
	  
	  *rans = *rans +	
	    2*(rW[(int) (bb+*nrW*gg)]*rW[(int) (aa+*nrW*dd)] +
	    rW[(int) (aa+*nrW*gg)]*rW[(int) (bb+*nrW*dd)]);
	}
      }
    } else { /* *ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;

	  *rans = *rans +	
	    2*(rW[(int) (aa+*nrW*gg)]*rW[(int) (bb+*nrW*gg)]);
	}
      ;
      }
    }
  } else { /* *ncA==1 */
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;
	  
	  *rans = *rans +	
	    2*(rW[(int) (aa+*nrW*gg)]* rW[(int) (aa+*nrW*dd)]);
	}
      }
    } else { /* *ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  
	  *rans = *rans +	
	    rW[(int) (aa+*nrW*gg)]*rW[(int) (aa+*nrW*gg)];
	}
      }
    }
  }


}


void trAWBVprim(double *rA, int *nrA, int *ncA,
		double *rW, int *nrW, int *ncW,
		double *rB, int *nrB, int *ncB, 
		double *rV, int *nrV, int *ncV, 
		double *rans)
{

  int i,j;
  int aa,bb,gg,dd;
  *rans=0;

  if (*ncA==2){
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;
	  //Rprintf(" %i %i %i %i\n", aa, bb, gg, dd);
	  *rans = *rans +	
	    rW[(int) bb+*nrV*gg] * rV[(int) dd+*nrW*aa] +
	    rW[(int) bb+*nrV*dd] * rV[(int) gg+*nrW*aa] +
	    rW[(int) aa+*nrV*gg] * rV[(int) dd+*nrW*bb] +
	    rW[(int) aa+*nrV*dd] * rV[(int) gg+*nrW*bb]; 
	}
      }
    } else { /* *ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	bb = (int)  rA[i + *nrA*1]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int) rB[j + *nrB*0]-1;	  
	  *rans = *rans +	
	    rV[(int) gg+*nrW*aa] * rW[(int) bb + *nrV*gg] +
	    rV[(int) gg+*nrW*bb] * rW[(int) aa + *nrV*gg];
	}
      }
    }
  } else { /* *ncA==1 */
    if (*ncB==2){
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;
	  dd = (int)  rB[j + *nrB*1]-1;  
	  *rans = *rans +	
	    rW[(int) aa+*nrV*gg] * rV[(int) dd+*nrW*aa] +
	    rW[(int) aa+*nrV*dd] * rV[(int) gg+*nrW*aa];
	}
      }
    } else { /* *ncB==1 */
      for (i=0; i<*nrA; i++){
	aa = (int)  rA[i + *nrA*0]-1;
	for (j=0; j<*nrB; j++){
	  gg = (int)  rB[j + *nrB*0]-1;	  
	  *rans = *rans +	
	    rW[(int) (aa+*nrV*gg)]*rV[(int) (gg+*nrW*aa)];
	}
      }
    }
  }
}

/* End of primitive functions */


SEXP trAW(SEXP A, SEXP W)
{
  int nrA, ncA, nrW, ncW;
  double *rans, *rA, *rW;
  SEXP Adims, Wdims, ans;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);

  PROTECT(A = AS_NUMERIC(A));
  PROTECT(W = AS_NUMERIC(W));

  rA = REAL(A);
  rW = REAL(W);

  nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1];
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 1000.0;

  trAWprim(rA, &nrA, &ncA,
	   rW, &nrW, &ncW, rans);

  //Rprintf("trAWB: %f\n",rans);
  UNPROTECT(3);
  return(ans);
}



SEXP trAWB(SEXP A, SEXP W, SEXP B)
{
  int nrA, ncA, nrW, ncW, nrB, ncB;
  double *rans, *rA, *rW, *rB;
  SEXP Adims, Wdims, Bdims, ans;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);
  Bdims = getAttrib(B, R_DimSymbol);


  PROTECT(A = AS_NUMERIC(A));
  PROTECT(W = AS_NUMERIC(W));
  PROTECT(B = AS_NUMERIC(B));

  rA = REAL(A);
  rW = REAL(W);
  rB = REAL(B);

  nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1];
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
  nrB = INTEGER(Bdims)[0];  ncB = INTEGER(Bdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 1000.0;

  trAWBprim(rA, &nrA, &ncA,
	    rW, &nrW, &ncW,
	    rB, &nrB, &ncB, rans);

  //Rprintf("trAWB: %f\n",rans);
  UNPROTECT(4);
  return(ans);
}


SEXP trAWBlist(SEXP Alist, SEXP W, SEXP Blist, SEXP mode)
{
  R_len_t nA = length(Alist), nB=length(Blist);
  SEXP ans;
  SEXP Aitem, Aidims, Bitem, Bidims, Wdims;
  int nrA, ncA, nrB, ncB, nrW, ncW, ii, jj;
  int Aii, Bjj, Astart, Aend, idx;
  double  *modep;

  double *rA, *rB, *rW, *ansp, rans;
  //if(!isNewList(Alist)) error("'list’ must be a list");

  Wdims = getAttrib(W, R_DimSymbol);
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
  rW  = REAL(W);

  PROTECT(ans = allocVector(REALSXP,nA*nB));
  ansp = REAL(ans);

  PROTECT(mode = coerceVector(mode, REALSXP)) ;
  modep = REAL(mode);

  idx = 0;

  for (Bjj=0; Bjj<nB; Bjj++){
    PROTECT(Bitem  = AS_NUMERIC(VECTOR_ELT(Blist, Bjj)));
    PROTECT(Bidims = getAttrib(Bitem, R_DimSymbol));
    
    nrB = INTEGER(Bidims)[0];  ncB = INTEGER(Bidims)[1];
    
    rB = REAL(Bitem);
    // printmatd(rB, &nrB, &ncB);

    if (*modep==0){
      Astart=0;      Aend  =nA;
    } else {
      Astart=Bjj;    Aend  =nA;
    }

    for (Aii=Astart; Aii<Aend; Aii++){
      
      PROTECT(Aitem  = AS_NUMERIC(VECTOR_ELT(Alist, Aii)));
      PROTECT(Aidims = getAttrib(Aitem, R_DimSymbol));
      
      nrA = INTEGER(Aidims)[0]; ncA = INTEGER(Aidims)[1];
      
      rA = REAL(Aitem);
      // printmatd(rA, &nrA, &ncA);
      
      trAWBprim(rA, &nrA, &ncA,
		rW, &nrW, &ncW,
		rB, &nrB, &ncB, &rans);

      //Rprintf("A: %i B: %i %f \n", Aii, Bjj, rans);
      ansp[idx++] = rans;
      UNPROTECT(2);
    }
    UNPROTECT(2);

  }
  
  //setAttrib(ans, R_NamesSymbol, getAttrib(Alist, R_NamesSymbol));
  UNPROTECT(2);
  return(ans);
}






SEXP trAWBW(SEXP A, SEXP W, SEXP B)
{
  int nrA, ncA, nrW, ncW, nrB, ncB;
  /*   int k1, k2; */
  int i, j;
  double *rans, *rA, *rW, *rB;
  SEXP Adims, Wdims, Bdims, ans;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);
  Bdims = getAttrib(B, R_DimSymbol);

  PROTECT(A = AS_NUMERIC(A));
  PROTECT(W = AS_NUMERIC(W));
  PROTECT(B = AS_NUMERIC(B));

  rA = REAL(A);
  rW = REAL(W);
  rB = REAL(B);

  nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1];
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
  nrB = INTEGER(Bdims)[0];  ncB = INTEGER(Bdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 0.0;

  trAWBWprim(rA, &nrA, &ncA,
	     rW, &nrW, &ncW,
	     rB, &nrB, &ncB, rans);

  //Rprintf("trAWBW: %f\n",rans);

  UNPROTECT(4);
  return(ans);
}




SEXP trAWBWlist(SEXP Alist, SEXP W, SEXP Blist, SEXP mode)
{
  R_len_t nA = length(Alist), nB=length(Blist);
  SEXP ans;
  SEXP Aitem, Aidims, Bitem, Bidims, Wdims;
  int nrA, ncA, nrB, ncB, nrW, ncW, ii, jj;
  int Aii, Bjj, Astart, Aend, idx;
  double  *modep;

  double *rA, *rB, *rW, *ansp, rans;
  //if(!isNewList(Alist)) error("'list’ must be a list");

  Wdims = getAttrib(W, R_DimSymbol);
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
  rW  = REAL(W);

  PROTECT(ans = allocVector(REALSXP,nA*nB));
  ansp = REAL(ans);

  PROTECT(mode = coerceVector(mode, REALSXP)) ;
  modep = REAL(mode);

  idx = 0;

  for (Bjj=0; Bjj<nB; Bjj++){
    PROTECT(Bitem  = AS_NUMERIC(VECTOR_ELT(Blist, Bjj)));
    PROTECT(Bidims = getAttrib(Bitem, R_DimSymbol));
    if (length(Bidims) < 2) error("Bad Bidims");
    
    nrB = INTEGER(Bidims)[0];  
    ncB = INTEGER(Bidims)[1];
    
    //Rprintf("B: %i %i %i\n", Bjj, nrB, ncB);  
    rB = REAL(Bitem);
    // printmatd(rB, &nrB, &ncB);

    if (*modep==0){
      Astart=0;      Aend  =nA;
    } else {
      Astart=Bjj;    Aend  =nA;
    }

    for (Aii=Astart; Aii<Aend; Aii++){
  
      PROTECT(Aitem  = AS_NUMERIC(VECTOR_ELT(Alist, Aii)));
      PROTECT(Aidims = getAttrib(Aitem, R_DimSymbol));
      
      nrA = INTEGER(Aidims)[0];  
      ncA = INTEGER(Aidims)[1];
      
      // Rprintf("A: %i %i %i\n", Aii, nrA, ncA);  
      rA = REAL(Aitem);
      
      trAWBWprim(rA, &nrA, &ncA,
		rW, &nrW, &ncW,
		rB, &nrB, &ncB, &rans);
      //Rprintf("%i %i %f\n", Aii, Bjj, rans);
      //ansp[Aii+nA*Bjj] = rans;
      //Rprintf("A: %i B: %i %f \n", Aii, Bjj, rans);
      ansp[idx++] = rans;
      UNPROTECT(2);
    }
    UNPROTECT(2);

  }
  
  //setAttrib(ans, R_NamesSymbol, getAttrib(Alist, R_NamesSymbol));
  UNPROTECT(2);
  return(ans);
}



SEXP trAWBV(SEXP A, SEXP W, SEXP B, SEXP V)
{
  int nrA, ncA, nrW, ncW, nrB, ncB, nrV, ncV;
  double *rans, *rA, *rW, *rB, *rV;
  SEXP Adims, Wdims, Bdims, Vdims, ans;

  Adims = getAttrib(A, R_DimSymbol);
  Wdims = getAttrib(W, R_DimSymbol);
  Bdims = getAttrib(B, R_DimSymbol);
  Vdims = getAttrib(V, R_DimSymbol);

  PROTECT(A = AS_NUMERIC(A));
  PROTECT(W = AS_NUMERIC(W));
  PROTECT(B = AS_NUMERIC(B));
  PROTECT(V = AS_NUMERIC(V));

  rA = REAL(A);
  rW = REAL(W);
  rB = REAL(B);
  rV = REAL(V);

  nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1];
  nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
  nrB = INTEGER(Bdims)[0];  ncB = INTEGER(Bdims)[1];
  nrV = INTEGER(Vdims)[0];  ncV = INTEGER(Vdims)[1];

  PROTECT(ans =allocVector(REALSXP,1));
  rans      = REAL(ans);
  *rans = 1000.0;

  trAWBVprim(rA, &nrA, &ncA,
	     rW, &nrW, &ncW,
	     rB, &nrB, &ncB, 
	     rV, &nrV, &ncV, rans);

  //Rprintf("trAWB: %f\n",rans);
  UNPROTECT(5);
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





/* /\*  Not used anymore *\/ */

/* SEXP trawb(SEXP A, SEXP W, SEXP B) */
/* { */
/*   int nrA, ncA, nrW, ncW, nrB, ncB; */
/*   /\*   int k1, k2; *\/ */
/*   int i, j; */
/*   double *rans, *rA=REAL(A), *rW=REAL(W), *rB=REAL(B); */
/*   SEXP Adims, Wdims, Bdims, ans; */

/*   /\* int a,b,g,d; *\/ */
/*   int aa,bb,gg,dd; */

/*   Adims = getAttrib(A, R_DimSymbol); */
/*   Wdims = getAttrib(W, R_DimSymbol); */
/*   Bdims = getAttrib(B, R_DimSymbol); */

/*   nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1]; */
/*   nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1]; */
/*   nrB = INTEGER(Bdims)[0];  ncB = INTEGER(Bdims)[1]; */

/*   PROTECT(ans =allocVector(REALSXP,1)); */
/*   rans      = REAL(ans); */
/*   *rans = 0.0; */

/*   if (ncA==2){ */
/*     if (ncB==2){ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	bb = (int)  rA[i + nrA*1]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
/* 	  dd = (int)  rB[j + nrB*1]-1; */
	  
/* 	  *rans = *rans + */
/* 	    rW[(int) (bb+nrW*gg)]*(aa==dd) + */
/* 	    rW[(int) (aa+nrW*gg)]*(bb==dd) + */
/* 	    rW[(int) (bb+nrW*dd)]*(aa==gg) + */
/* 	    rW[(int) (aa+nrW*dd)]*(bb==gg) ; */
/* 	} */
/*       } */
/*     } else { /\* ncB==1 *\/ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	bb = (int)  rA[i + nrA*1]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */

/* 	  *rans = *rans + */
/* 	    rW[(int) (gg+nrW*aa)]*(gg==bb) + */
/* 	    rW[(int) (gg+nrW*bb)]*(gg==aa); */
/* 	} */
/*       ; */
/*       } */
/*     } */
/*   } else { /\* ncA==1 *\/ */
/*     if (ncB==2){ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
/* 	  dd = (int)  rB[j + nrB*1]-1; */
	  
/* 	  *rans = *rans + */
/* 	    rW[(int) (aa+nrW*gg)]*(aa==dd) + */
/* 	    rW[(int) (aa+nrW*dd)]*(aa==gg); */
/* 	} */
/*       } */
/*     } else { /\* ncB==1 *\/ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
	  
/* 	  *rans = *rans + */
/* 	    rW[(int) (aa+nrW*aa)]*(aa==gg); */
/* 	} */
/*       } */

/*     } */
/*   } */
  
/*   UNPROTECT(1); */
/*   return(ans); */
/* } */





/* Not used anymore */

/* SEXP trawbw(SEXP A, SEXP W, SEXP B) */
/* { */
/*   int nrA, ncA, nrW, ncW, nrB, ncB; */
/*   /\*   int k1, k2; *\/ */
/*   int i, j; */
/*   double *rans, *rA=REAL(A), *rW=REAL(W), *rB=REAL(B); */
/*   SEXP Adims, Wdims, Bdims, ans; */

/*   /\*   int a,b,g,d; *\/ */
/*   int aa,bb,gg,dd; */

/*   Adims = getAttrib(A, R_DimSymbol); */
/*   Wdims = getAttrib(W, R_DimSymbol); */
/*   Bdims = getAttrib(B, R_DimSymbol); */

/*   nrA = INTEGER(Adims)[0]; */
/*   ncA = INTEGER(Adims)[1]; */
  
/*   nrW = INTEGER(Wdims)[0]; */
/*   ncW = INTEGER(Wdims)[1]; */

/*   nrB = INTEGER(Bdims)[0]; */
/*   ncB = INTEGER(Bdims)[1]; */

/*   PROTECT(ans =allocVector(REALSXP,1)); */
/*   rans      = REAL(ans); */
/*   *rans = 0.0; */

/*   if (ncA==2){ */
/*     if (ncB==2){ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	bb = (int)  rA[i + nrA*1]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
/* 	  dd = (int)  rB[j + nrB*1]-1; */
	  
/* 	  *rans = *rans +	 */
/* 	    2*(rW[(int) (bb+nrW*gg)]*rW[(int) (aa+nrW*dd)] + */
/* 	    rW[(int) (aa+nrW*gg)]*rW[(int) (bb+nrW*dd)]); */
/* 	} */
/*       } */
/*     } else { /\* ncB==1 *\/ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	bb = (int)  rA[i + nrA*1]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */

/* 	  *rans = *rans +	 */
/* 	    2*(rW[(int) (aa+nrW*gg)]*rW[(int) (bb+nrW*gg)]); */
/* 	} */
/*       ; */
/*       } */
/*     } */
/*   } else { /\* ncA==1 *\/ */
/*     if (ncB==2){ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
/* 	  dd = (int)  rB[j + nrB*1]-1; */
	  
/* 	  *rans = *rans +	 */
/* 	    2*(rW[(int) (aa+nrW*gg)]* rW[(int) (aa+nrW*dd)]); */
/* 	} */
/*       } */
/*     } else { /\* ncB==1 *\/ */
/*       for (i=0; i<nrA; i++){ */
/* 	aa = (int)  rA[i + nrA*0]-1; */
/* 	for (j=0; j<nrB; j++){ */
/* 	  gg = (int)  rB[j + nrB*0]-1; */
	  
/* 	  *rans = *rans +	 */
/* 	    rW[(int) (aa+nrW*gg)]*rW[(int) (aa+nrW*gg)]; */
/* 	} */
/*       } */
/*     } */
/*   } */
  
/*   UNPROTECT(1); */
/*   return(ans); */
/* } */
