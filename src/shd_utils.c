#include <Rdefines.h>
#include <R.h>
#include <stdlib.h>
#include <R_ext/Lapack.h> 
#include "shd_print.h";

/* *** Forward declarations *** */

/* Matrix operations */

void schursubt(double *K, int *nvar, int *gset, int *gsetlen, double *ans);

void schursubt2(double *K, int *nvar, 
		int *gset, int *gsetlen, 
		int *cset, int *csetlen, 
 		double *TMP, double *ans);

void shddet(double *Ain, int *nrA, double *ans);

void shdinverse(double *A, int *nrA, double *tolin);

void shdmatprod(double *x, int *nrx, int *ncx,
		double *y, int *nry, int *ncy, double *z);

void shdsolve(double *A, int *nrA, double *B, int *ncB, double *tolin);

void shdsubmatrix(double *A, int *nrA, int *ncA, 
		  int *gen1, int *gen1len, int *gen2, int *gen2len, double *ans);

void shdtraceAB(double *A, int *nrA, int *ncA,
		double *B, int *nrB, int *ncB, double *ans);


/* Utilities */

void complement(int *set, int *setlen, int *nvar, int *ans);

void shdunique(double *gen, int *ngen, int *nunique, int *nvar, double *ans); 



/* Find K(gc)inv(K(cc))K(cc) */

void schursubt(double *K, int *nvar, int *gset, int *gsetlen, double *ans)
{
  int *cset, ng, nc;
  double *Awork, *Bwork, *Cwork, tol;

  ng = *gsetlen;
  nc = *nvar - *gsetlen;
  cset    = (int *) R_alloc(*nvar, sizeof(int));
  complement(gset, gsetlen, nvar, cset);
  //printveci(cset, nvar);

  Awork = (double*) malloc (nc*ng * sizeof(double)); 
  Bwork = (double*) malloc (nc*nc * sizeof(double)); 
  Cwork = (double*) malloc (ng*nc * sizeof(double)); 

      // K(cg) -> Awork
      //
      shdsubmatrix(K, nvar, nvar, cset, &nc, gset, &ng, Awork);
      //Rprintf("K(cg):\n"); printmatd(Awork, &nc, &ng);
      
      // inv(K(cc)) -> Bwork
      //
      shdsubmatrix(K, nvar, nvar, cset, &nc, cset, &nc, Bwork);
      //Rprintf("K(cc):\n"); printmatd(Bwork, &nc, &nc);

      // inv(K(cc))K(cg) -> Awork
      shdsolve(Bwork, &nc, Awork, &ng, &tol); 
      //Rprintf("inv(K(cc))K(cg):\n"); printmatd(Awork, &nc, &ng);

      // K(gc) -> Cwork
      //
      shdsubmatrix(K, nvar, nvar, gset, &ng, cset, &nc, Cwork);
      //Rprintf("K(gc):\n"); printmatd(Cwork, &ng, &nc);

      // K(gc)inv(K(cc))K(cg) -> ans
      //
      shdmatprod(Cwork, &ng, &nc, Awork, &nc, &ng, ans);    
      //Rprintf("K(gc)inv(K(cc))K(cg):\n"); printmatd(ans, &ng, &ng);
      
      free(Awork); free(Cwork);  free(Bwork);
}





/*
  Input: 
  K: *nvar**nvar matrix
  gset (g): A subset of 1:*nvar.
  cset (c): Complement of g
  Output:
  Partitions K as Kgg, Kcg, Kgc, Kcc; 
*/

void schurpart(double *K, int *nvar, int *gset, int *gsetlen, int *cset, int *csetlen, double *ans)
{
  int    ng=*gsetlen, nc=*csetlen;
  double *Kgg, *Kcg, *Kgc, *Kcc;

  Kgg = ans;
  Kcg = ans + ng*ng;
  Kgc = Kcg + nc*ng;
  Kcc = Kgc + ng*nc;

  shdsubmatrix(K, nvar, nvar, gset, &ng, gset, &ng, Kgg);
  //Rprintf("Kgg:\n"); printmatd(Kgg, &ng, &ng); 
  shdsubmatrix(K, nvar, nvar, cset, &nc, gset, &ng, Kcg);
  //Rprintf("Kcg:\n"); printmatd(Kcg, &nc, &ng); 
  shdsubmatrix(K, nvar, nvar, gset, &ng, cset, &nc, Kgc);
  //Rprintf("Kgc:\n"); printmatd(Kgc, &ng, &nc); 
  shdsubmatrix(K, nvar, nvar, cset, &nc, cset, &nc, Kcc);
  //Rprintf("Kcc:\n"); printmatd(Kcc, &nc, &nc); 
  //Rprintf("ans:\n"); printmatd(ans, nvar, nvar); 
}



/*
  TMP: Vector of same length (at least) as K
*/
void schursubt2(double *K, int *nvar, 
		int *gset, int *gsetlen, 
		int *cset, int *csetlen, 
 		double *TMP, double *ans)
{
  int    ng, nc;
  double tol;

  ng      = *gsetlen;
  nc      = *csetlen; // *nvar - *gsetlen;

  schurpart(K, nvar, gset, &ng, cset, &nc, TMP);

  double *Kgg, *Kcg, *Kgc, *Kcc;

  Kgg = TMP;
  Kcg = TMP + ng*ng;
  Kgc = Kcg + nc*ng;
  Kcc = Kgc + ng*nc;

  shdsolve(Kcc, &nc, Kcg, &ng, &tol); 
  //Rprintf("inv(K(cc))K(cg):\n"); printmatd(Kcg, &nc, &ng);

  shdmatprod(Kgc, &ng, &nc, Kcg, &nc, &ng, ans);   
  //Rprintf("K(gc)inv(K(cc))K(cg):\n"); printmatd(ans, &ng, &ng);
}
















/* Matrix operations - Implementation */

void shddet(double *Ain, int *nrA, double *ans)
{
  int i, n, info, *jpvt, sign, useLog=1;
  double modulus = 0.0; /* -Wall */
  int nrA2;
  double *Awork;
  
  nrA2 = *nrA * *nrA;
  
  Awork = (double *) R_alloc(*nrA * *nrA, sizeof(double));  
  Memcpy(Awork, Ain, (size_t) (nrA2));     
  
  n = *nrA;
  //jpvt = (int *) R_alloc(*nrA, sizeof(int));
  
  jpvt = (int *) malloc(*nrA *  sizeof(int));
  F77_CALL(dgetrf)(&n, &n, Awork, &n, jpvt, &info);
  
  //printmatd(Awork, &n, &n);
  
  sign = 1;
  if (info < 0)
    Rprintf("error code %d from Lapack routine '%s'", info, "dgetrf");
  else if (info > 0) { /* Singular matrix:  U[i,i] (i := info) is 0 */
    /*warning("Lapack dgetrf(): singular matrix: U[%d,%d]=0", info,info);*/
    modulus = (useLog ? R_NegInf : 0.);
  }
  else {
    for (i = 0; i < n; i++) if (jpvt[i] != (i + 1))
      sign = -sign;
    if (useLog) {
      modulus = 0.0;
      for (i = 0; i < n; i++) {
	double dii = Awork[i*(n + 1)]; /* ith diagonal element */
	modulus += log(dii < 0 ? -dii : dii);
	if (dii < 0) sign = -sign;
      }
    } else {
      modulus = 1.0;
      for (i = 0; i < n; i++)
	modulus *= Awork[i*(n + 1)];
      if (modulus < 0) {
	modulus = -modulus;
	sign = -sign;
      }
    }
  }
  
  //Rprintf("%i %f %f \n", sign, modulus, ((float) sign) * exp(modulus));
  *ans = sign * exp(modulus);
  free(jpvt);
}



/* shdinverse
   ==========
   Finds the inverse of A; returns the answer in A
*/

void shdinverse(double *A, int *nrA, double *tolin)
{
  int ii, jj, info, *ipiv, n=*nrA, n2=*nrA**nrA;
  double *B;

  /*   Rprintf("shdinverse (entry):\n");  */
  /*   Rprintf("nrA: %i\n", *nrA); */
  /*   printmatd(A, nrA, nrA);  */
  B     = (double *) R_alloc(n2, sizeof(double)); //Rprintf("B ok...\n"); 
  ipiv  = (int *)    R_alloc(n, sizeof(int));     //Rprintf("ipiv ok...\n"); 

  //B     = (double *) malloc(n2 * sizeof(double)); //Rprintf("B ok...\n");
  //ipiv  = (int *)    malloc(n *sizeof(int)); //Rprintf("ipiv ok...\n");
  
  for (ii=0; ii<*nrA; ii++){
    for (jj=0; jj<*nrA; jj++){
      if (ii==jj){
	B[ii+*nrA*jj] = 1;
      } else {
	B[ii+*nrA*jj] = 0;
      }
    }
  }
  //Rprintf("diag:\n"); printmatd(B, nrA, nrA);  

  F77_CALL(dgesv)(nrA, nrA, A, nrA, ipiv, B, nrA, &info);

  Memcpy(A, B, (size_t) n2);   
  //Rprintf("shdinverse (exit):\n"); printmatd(A, nrA, nrA);
  //free(ipiv); free(B);
}


/* 
   Calculates product of matrices x and y; returns result in z 
*/

void shdmatprod(double *x, int *nrx, int *ncx,
		double *y, int *nry, int *ncy, double *z)
{
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, nrx, ncy, ncx, &one,
		  x, nrx, y, nry, &zero, z, nrx);
}


/* shdsolve
   ========
   Solves the system Ax=B; returns the answer in B
*/

void shdsolve(double *A, int *nrA, double *B, int *ncB, double *tolin)
{
  int info;
  //int *ipiv  = (int *) R_alloc((int) *nrA, sizeof(int));    
  int *ipiv  = (int *) malloc (*nrA * sizeof(int)); 
  F77_CALL(dgesv)(nrA, ncB, A, nrA, ipiv, B, nrA, &info);
  free(ipiv);
}



/* shdsubmatrix:
   ==========
   A: On entry: A is the matrix from which a submatrix (defined by gen1 and gen2)  
   must be found 
   A: On exit: A contains the submatrix defined by gen1, gen2 in the first gen1len*gen2len 
   entries 
*/

void shdsubmatrix(double *A, int *nrA, int *ncA, 
	       int *gen1, int *gen1len, int *gen2, int *gen2len, double *ans)
{
  int ii, jj, kk;  
  //Rprintf("In shdsubmatrix:\n");
  //printveci(gen1, gen1len);  printveci(gen2, gen2len); printmatd(A, nrA, ncA);

  kk = 0;
  for (jj=0; jj<*gen2len; jj++){
    for (ii=0; ii<*gen1len; ii++){
      // Rprintf("%i %i %i %i \n", ii, gen1[ii], jj, gen2[jj]);
      ans[kk++] = A[gen1[ii] + *nrA * gen2[jj]];
    }
  }
} 

void shdtraceAB(double *A, int *nrA, int *ncA,
		double *B, int *nrB, int *ncB, double *ans)
{
  int ii;
  double x=0;

  for (ii=0; ii<*nrA**nrA; ii++){
    x=x+A[ii]*B[ii];
  }
  
  //Rprintf("tr: %f\n", x);
  *ans=x;
}




/* Utilities - Implementation */

/* 
   Find the complement of set in 0:nvar-1 
   E.g. set=1,2,4,5; nvar=10; returns  0 3 6 7 8 9 1 2 4 5 
   so complement is at beginning, set is at the end and both are sorted. 
*/
void complement(int *set, int *setlen, int *nvar, int *ans)
{
  int ii,jj=0, kk=0, *set2;
  set2 = (int*) malloc (*setlen * sizeof(double)); 

  for (ii=0; ii<*nvar; ii++)
    ans[ii] = ii;
  //printveci(ans, nvar);

  for (ii=0; ii<*setlen; ii++)
    ans[set[ii]] = -1;

  //printveci(ans, nvar);

  for (ii=0; ii<*nvar; ii++){
    if (ans[ii]>=0){
      ans[kk++] = ans[ii];
    } else {
      set2[jj++] = ii;
    }
  }
  
  jj = 0;
  for (ii=kk; ii<*nvar; ii++){
    ans[ii] = set2[jj++];
  }
  free(set2);
  //printveci(ans, nvar);
}



/* 
   On entry:
   *gen :
   vector with *ngen elements, e.g. gen=(2,5,2,45,23,45), 
   *ngen :
   length of *gen, e.g. *ngen=6.
   *nunique :
   placeholder for answer info
   *nvar :
   the maximal possible element size in gen, e.g. 100; 
   *ans : 
   vector with at least *nvar elements. 
   
   On exit:
   *ans contains (2,5,23,45) and *nunique is 4
*/

int doublecmp(const void *v1, const void *v2){
  return (*(double *)v1 - *(double *)v2);
}

void shdunique(double *gen, int *ngen, int *nunique, int *nvar, double *ans)
{
  int ii, jj, found, kk=0;
  if (*ngen==1){
    kk=*nunique = 1;
    ans[0] = gen[0];
    //Rprintf("%f %f\n", ans[0], gen[0]); printvecd(ans, nvar);
  } else { 
    Memcpy(ans, gen, (size_t) *ngen);              // printvecd(ans, ngen);
    qsort(ans, *ngen, sizeof(double), doublecmp);  // printvecd(ans, ngen);
    kk = 1;
    for (ii=1; ii<*ngen;ii++){
      if (ans[ii]!=ans[ii-1]){
	ans[kk++] = ans[ii];
      }
    }
    *nunique = kk;
    //for (ii=kk; ii<*nvar;ii++)
    //  ans[ii] = -1;
  }
  
  for (ii=0; ii<*nvar; ii++){
    found=0;
    for (jj=0; jj<*nunique; jj++){
      if (ans[jj]==ii){
	found=1;
	break;
      }
    }
    if (found==0){
      //Rprintf(" %i\n", ii);
      ans[kk++] = ii;
      
    }
  }

  //Rprintf(" %i\n", *nunique);  
}







/* void printmatf(float *Avec, int *nr, int *nc){ */
/*   int ii, jj; */
  
/*   for (ii=0; ii<*nr; ii++){ */
/*     for (jj=0; jj<*nc; jj++){ */
/*       Rprintf(" %5.2f ", Avec[ii + *nr * jj]); */
/*     } */
/*     Rprintf("\n"); */
/*   } */
/*   Rprintf(" ---------------------------\n"); */
/* } */



    //Adims = INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP));
    //n = Adims[0];
    //if (Adims[1] != n)
    //error(_("'a' must be a square matrix"));
    //jpvt = (int *) R_alloc(n, sizeof(int));


    //if (!(isMatrix(Ain) && isReal(Ain)))
    //error(_("'a' must be a numeric matrix"));
    //useLog = asLogical(logarithm);
    //if (useLog == NA_LOGICAL) error(_("argument 'logarithm' must be logical"));
    //PROTECT(A = duplicate(Ain));




/* void submatrix(double *A, int *nrA, int *ncA,  */
/* 	       int *gen1, int *gen1len, int *gen2, int *gen2len) */



/*   // printmatd(S, nvar, nvar); printmati(gmat, ngen, nvar); */

/*   Saa = (double *) R_alloc(*nvar * *nvar, sizeof(double)); */
/*   Inn = (double *) R_alloc(*nvar * *nvar, sizeof(double));  */
/*   Iaa = (double *) R_alloc(*nvar * *nvar, sizeof(double));  */


/*   for (ii=0; ii<*nvar; ii++){ */
/*     for (jj=0; jj<*nvar; jj++){ */
/*       if (ii==jj) */
/* 	Inn[ii + *nvar *jj] = 1; */
/*       else */
/* 	Inn[ii + *nvar *jj] = 0; */
/*     } */
/*   } */
/*   printmatd(Inn, nvar, nvar); */

  
/*   currg = 0;   */
  
/*   kk = 0; */
/*   for (ii=gs; ii<ge; ii++){  */
/*     ii2 = gmat[currg + *ngen * ii]-1; */
/*     for (jj=gs; jj<ge; jj++){ */
/*       jj2 = gmat[currg + *ngen * jj]-1; */
/*       // Rprintf("%i %i %f \n", ii2, jj2, S[ii2 + *nvar * jj2]); */
/*       Saa[kk] = S[ii2 + *nvar * jj2]; */
/*       Iaa[kk] = Inn[ii2 + *nvar * jj2]; */
/*       kk++; */
/*     } */
/*   } */
  
/*   // Saa er (ge x ge) matrix; Find den inverse */
/*   printmatd(Saa, &ge, &ge);  printmatd(Iaa, &ge, &ge); */

/*   ipiv  = (int *) R_alloc(ge, sizeof(int)); */

/*   F77_CALL(dgesv)(nvar, &ge, Saa, nvar, ipiv, Iaa, nvar, &info); */
/*   if (info == 0) */
/*     Rprintf("ok...\n"); */


/*   printmatd(Saa, &ge, &ge);   */
/*   printmatd(Iaa, &ge, &ge); */



/* 
   gc and gclen defines {a1, a2, ..., ap} 
   n (derived) and {ai} defines complements {b1,...bp}

   We need A[a,a]; inv(A[a,a]); A[a,b]


*/









/*   gs = 0; */
/*   for (kk=0; kk< (int) ng; kk++){ */
/*     for (ii=gs; ii< (int) gclenp[kk]; ii++){ */

/*     } */
/*   } */


/*   ge = 2; */
/*   Aaan = ge-gs+1; */
/*   Aaax = (double *) R_alloc(Aaan * Aaan, sizeof(double)); */
  
/*   for (ii=gs; ii<=ge; ii++){ */
/*     Rprintf(" %i \n", (int) gcp[ii]); */
/*   } */

/*   kk=0; */
/*   for (ii=gs; ii<=ge; ii++){  */
/*     for (jj=gs; jj<=ge; jj++){ */
/*       Aaax[kk++] = Ax[(int)(gcp[jj]+n*gcp[ii])]; */
/*     } */
/*   } */

/*   printmat(Aaax, &Aaan, &Aaan); */






      /* Rprintf("%i %i %i %f\n", (int)gcp[ii], (int)gcp[jj], (int)(gcp[ii]+n*gcp[jj]), 
	 avals[(int)(gcp[ii]+n*gcp[jj])] ); */



/*   anorm = F77_CALL(dlange)("1", &n, &n, REAL(A), &n, (double*) NULL); */
/*   work = (double *) R_alloc(4*n, sizeof(double)); */
/*   F77_CALL(dgecon)("1", &n, avals, &n, &anorm, &rcond, work, ipiv, &info); */
/*   if (rcond < tol) */
/*     Rprintf("system is computationally singular: reciprocal condition number = %g", */
/* 	    rcond); */


/*   if(n == 0) error("'a' is 0-diml"); */
/*   if(p == 0) error("no right-hand side in 'b'"); */
/*   if(Adims[1] != n) */
/*     error(_("'a' (%d x %d) must be square"), n, Adims[1]); */
/*   if(Bdims[0] != n) */
/*     error(_("'b' (%d x %d) must be compatible with 'a' (%d x %d)"), */
/* 	  Bdims[0], p, n, n); */


/*   if (info < 0) */
/*     error(_("argument %d of Lapack routine %s had invalid value"), */
/* 	  -info, "dgesv"); */
/*   if (info > 0) */
/*     error(_("Lapack routine dgesv: system is exactly singular")); */



/*     for (ii=0; ii<ng; ii++){ */
/*       for (jj=0; jj<ng; jj++){ */
/* 	sum = 0;     */
/* 	for (kk=0; kk<nc; kk++){ */
/* 	  sum = sum + Awork[ii+nc*kk] * Cwork[kk+nc*jj]; */
/* 	} */
/* 	//Rprintf("%i %i sum %f, S %f\n", ii, jj, sum, Swork[ii+ng*jj]); */
/* 	Bwork[ii+ng*jj] = sum ; //+ Swork[ii+ng*jj]; */
/*       } */
/*     } */
    
/*     Rprintf("B (iteration):\n"); printmatd(Bwork, &ng, &ng); */


/*     for (ii=0; ii<nc; ii++){ */
/*       for (jj=0; jj<ng; jj++){ */
/* 	sum = 0; */
/* 	for (kk=0; kk<nc; kk++){ */
/* 	  sum = sum + Awork[ii + nc * kk] *Bwork[kk + nc * jj]; */
/* 	} */
/* 	Cwork[ii+nc*jj] = sum; */
/*       } */
/*     } */
/*     Rprintf("inv(K(cc))K(cg):\n"); printmatd(Cwork, &nc, &ng); */


    //shdinverse(K, &clen[currg], &tolin);
    //Rprintf("Its inverse:\n"); printmatd(K, &nc, &nc);

    // Rprintf("K (iteration):\n"); printmatd(K, nvar, nvar);

