#include <Rinternals.h>



/* void printArray(int *array, int nr, int nc) {  */
/*   int i, j; */
  
/*   Rprintf("JJJJJJJJJJJJJ\n"); */
/*   for(i=0 ; i<nr ; i++) { */
/*     for(j=0 ; j<nc ; j++) { */
/*       Rprintf("%2d ", array[i*nr+j]); */
/*     } */
/*     Rprintf("\n"); */
/*   } */
/* } */

SEXP eim(SEXP x, SEXP y)
{
  int nrx, ncx, nry, ncy, mode;
  int i, j, kk=0;
  double *rans, *rx=REAL(x), *ry=REAL(y);
  SEXP xdims, ydims, ans;
  xdims = getAttrib(x, R_DimSymbol);
  ydims = getAttrib(y, R_DimSymbol);
  nrx = INTEGER(xdims)[0];
  ncx = INTEGER(xdims)[1];
  nry = INTEGER(ydims)[0];
  ncy = INTEGER(ydims)[1];
  
  int xa[nrx][ncx];
  for (i=0; i<nrx; i++){
    for (j=0; j<ncx; j++){
      xa[i][j] = rx[i+nrx*j];
      /*Rprintf(" %i", xa[i][j]);*/
    }
    /* Rprintf("\n");*/
  }


  int ya[nry][ncy];
  for (i=0; i<nry; i++){
    for (j=0; j<ncy; j++){
      ya[i][j] = ry[i+nry*j];
      /*Rprintf(" %i", ya[i][j]);      */
    }
    /* Rprintf("\n");*/
  }


  /*   Initialize */
  int nrans = nrx*nry*2;  
  int ansa[nrans][4];
  /*
  for (i=0; i<nrans; i++){
    for (j=0; j<4; j++){
      ansa[i][j] = 0;
    }
  }
  */

  kk = 0; 
  if (ncx==2){
    for (i=0;i<nrx;i++){
      for (j=0;j<nry*2;j++){
	ansa[kk][0] = xa[i][0];
	ansa[kk][2] = xa[i][1];
	kk++;
      }
    }
  } else {
    for (i=0;i<nrx;i++){
      for (j=0;j<nry*2;j++){
	ansa[kk][0] = xa[i][0];
	ansa[kk][2] = xa[i][0];
	kk++;
      }
    } 
  }


/*   for (i=0; i<nrans; i++){ */
/*     for (j=0; j<4; j++){ */
/*       Rprintf(" %i ", ansa[i][j]);       */
/*     } */
/*     Rprintf("\n"); */
/*   } */
  
  kk = 0;
  if (ncy==2){
    for (i=0;i<nrx;i++){
      for (j=0;j<nry;j++){
	ansa[kk][1] = ya[j][1];
	ansa[kk][3] = ya[j][0];
	kk++;
	ansa[kk][1] = ya[j][0];
	ansa[kk][3] = ya[j][1];
	kk++;
      }
    }
  } else {
    for (i=0;i<nrx;i++){
      for (j=0;j<nry;j++){
	ansa[kk][1] = ya[j][0];
	ansa[kk][3] = ya[j][0];
	kk++;
      }
      for (j=0;j<nry;j++){
	ansa[kk][1] = ya[j][0];
	ansa[kk][3] = ya[j][0];
	kk++;
      }    
    }
  }
  
/*   for (i=0; i<nrans; i++){ */
/*     for (j=0; j<4; j++){ */
/*       Rprintf(" %i ", ansa[i][j]);       */
/*     } */
/*     Rprintf("\n"); */
/*   } */

  mode      = REALSXP;
  PROTECT(ans = allocMatrix(mode, nrans, 4));
  rans      = REAL(ans);

  for (i=0; i<nrans; i++){
    for (j=0; j<4; j++){
      rans[i+nrans*j] = ansa[i][j];
    }
  }

  PROTECT(ans = coerceVector(ans, INTSXP));
  UNPROTECT(2);
  return(ans);
}



/*   PROTECT(tr   = allocVector(mode,1)); */

/*   REAL(tr)[0] = 0; */
/*   matprod(REAL(x), nrx, ncx,  */
/*  	  REAL(y), nry, ncy, REAL(ans));  */
  
/*   for (i=0; i< nrx; i++){ */
/*     REAL(tr)[0] = REAL(tr)[0] + REAL(ans)[i*(nrx+1)]; */
/*   } */

/*   Rprintf("%f", REAL(tr)[0]); */


/* static void eim(double *x, int nrx, int ncx, */
/* 		double *y, int nry, int ncy,  */
/* 		double *z, int nrz, int ncz) */
/* { */
/* /\*   char *transa = "N", *transb = "N"; *\/ */
/* /\*   double one = 1.0, zero = 0.0; *\/ */
/* /\*   F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one, *\/ */
/* /\* 		    x, &nrx, y, &nry, &zero, z, &nrx); *\/ */

/* } */

