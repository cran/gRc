#include <RcppArmadillo.h>
// #include <Rdefines.h>
// #include <R.h>

using namespace Rcpp;
//[[Rcpp::interfaces(r,cpp)]]

//[[Rcpp::export]]
double trAW_(NumericMatrix rA, NumericMatrix rW)
{

  double ans=0;
  int nrA = rA.nrow(), ncA = rA.ncol();
  int nrW = rW.nrow(); //, ncW = rW.ncol();

  int ii, idx;

  if (ncA == 2){
    for (ii=0; ii<nrA; ii++){
      idx = (int) (rA[ii] - 1 + nrW * (rA[ii + nrA] - 1));
      //Rprintf("%i ", idx);
      ans = ans + rW[idx];
    }
    ans = 2 * ans;    
  } else {
    for (ii=0; ii < nrA; ii++){
      idx = (int) ((rA[ii] - 1) * (nrW + 1));
      //Rprintf("%i %f \n", idx, rA[ii]);
      ans = ans + rW[idx];
    }
  }
  return(ans);
}

//[[Rcpp::export]]
double trAWB_(NumericMatrix rA, NumericMatrix rW, NumericMatrix rB)
{

  double ans=0;
  int nrA = rA.nrow(), ncA = rA.ncol();
  int nrW = rW.nrow(); //, ncW = rW.ncol();
  int nrB = rB.nrow(), ncB = rB.ncol();
  
  int i,j;
  int aa,bb,gg,dd;

  
  //printmatd(rA, nrA, ncA);
  //printmatd(rB, nrB, ncB);
  
  if (ncA == 2){
    if (ncB == 2){
      for (i=0; i < nrA; i++){
	aa = (int)  rA[i] - 1;
	bb = (int)  rA[i + nrA] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int)  rB[j] - 1;
	  dd = (int)  rB[j + nrB] - 1;
	  
	  ans = ans +
	    rW[(bb + nrW * gg)] * (aa==dd) +
	    rW[(aa + nrW * gg)] * (bb==dd) +
	    rW[(bb + nrW * dd)] * (aa==gg) +
	    rW[(aa + nrW * dd)] * (bb==gg) ;
	}
      }
    } else { /* ncB==1 */
      for (i=0; i < nrA; i++){
	aa = (int)  rA[i] - 1;
	bb = (int)  rA[i + nrA] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int) rB[j] - 1;
	  
	  ans = ans +
	    rW[(gg + nrW * aa)] * (gg==bb) +
	    rW[(gg + nrW * bb)] * (gg==aa);
	}
	;
      }
    }
  } else { /* ncA==1 */
    if (ncB == 2){
      for (i=0; i < nrA; i++){
	aa = (int) rA[i] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int) rB[j] - 1;
	  dd = (int) rB[j + nrB] - 1;
	  
	  ans = ans +
	    rW[(aa + nrW * gg)] * (aa==dd) +
	    rW[(aa + nrW * dd)] * (aa==gg);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i < nrA; i++){
	aa = (int) rA[i] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int) rB[j] - 1;
	  
	  ans = ans +
	    rW[(aa + nrW * aa)] * (aa==gg);
	}
      }

    }
  }
  return(ans);
}


//[[Rcpp::export]]
double trAWBW_(NumericMatrix rA, NumericMatrix rW, NumericMatrix rB)
{

  double ans=0;
  int nrA = rA.nrow(), ncA = rA.ncol();
  int nrW = rW.nrow(); //, ncW = rW.ncol();
  int nrB = rB.nrow(), ncB = rB.ncol();

  int i,j;
  int aa,bb,gg,dd;

  if (ncA == 2){
    if (ncB == 2){
      for (i=0; i < nrA; i++){
	aa = (int)  rA[i] - 1;
	bb = (int)  rA[i + nrA] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int)  rB[j] - 1;
	  dd = (int)  rB[j + nrB] - 1;
	  
	  ans = ans +	
	    2 * (rW[(bb + nrW * gg)] * rW[(aa + nrW * dd)] +
		 rW[(aa + nrW * gg)] * rW[(bb + nrW * dd)]);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i <nrA; i++){
	aa = (int)  rA[i] - 1;
	bb = (int)  rA[i + nrA] - 1;
	for (j=0; j <nrB; j++){
	  gg = (int) rB[j] - 1;

	  ans = ans +	
	    2 * (rW[(aa + nrW * gg)] * rW[(bb + nrW * gg)]);
	}
      ;
      }
    }
  } else { /* ncA==1 */
    if (ncB == 2){
      for (i=0; i < nrA; i++){
	aa = (int)  rA[i] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int)  rB[j] - 1;
	  dd = (int)  rB[j + nrB] - 1;
	  
	  ans = ans +	
	    2*(rW[(aa + nrW * gg)] * rW[(aa + nrW * dd)]);
	}
      }
    } else { /* ncB==1 */
      for (i=0; i<nrA; i++){
	aa = (int) rA[i] - 1;
	for (j=0; j < nrB; j++){
	  gg = (int)  rB[j] - 1;
	  
	  ans = ans +	
	    rW[(aa + nrW * gg)] * rW[(aa + nrW * gg)];
	}
      }
    }
  }


  
  return(ans);
}


//[[Rcpp::export]]
double trAWBV_(NumericMatrix rA, NumericMatrix rW, NumericMatrix rB, NumericMatrix rV)
{

  double ans=0;
  int nrA = rA.nrow(), ncA = rA.ncol();
  int nrW = rW.nrow(); //, ncW = rW.ncol();
  int nrB = rB.nrow(), ncB = rB.ncol();
  int nrV = rV.nrow(); //, ncV = rV.ncol();
  
  int i,j;
  int aa,bb,gg,dd;

  if (ncA == 2){
    if (ncB == 2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i ]-1;
	bb = (int)  rA[i + nrA]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j]-1;
	  dd = (int)  rB[j + nrB]-1;
	  //Rprintf(" %i %i %i %i\n", aa, bb, gg, dd);
	  ans = ans +	
	    rW[bb + nrV * gg] * rV[dd + nrW * aa] +
	    rW[bb + nrV * dd] * rV[gg + nrW * aa] +
	    rW[aa + nrV * gg] * rV[dd + nrW * bb] +
	    rW[aa + nrV * dd] * rV[gg + nrW * bb]; 
	}
      }
    } else { /* ncB == 1 */
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i]-1;
	bb = (int)  rA[i + nrA]-1;
	for (j=0; j<nrB; j++){
	  gg = (int) rB[j]-1;	  
	  ans = ans +	
	    rV[gg + nrW * aa] * rW[bb + nrV * gg] +
	    rV[gg + nrW * bb] * rW[aa + nrV * gg];
	}
      }
    }
  } else { /* ncA == 1 */
    if (ncB == 2){
      for (i=0; i<nrA; i++){
	aa = (int)  rA[i]-1;
	for (j=0; j<nrB; j++){
	  gg = (int)  rB[j]-1;
	  dd = (int)  rB[j + nrB]-1;  
	  ans = ans +	
	    rW[aa + nrV * gg] * rV[dd + nrW * aa] +
	    rW[aa + nrV * dd] * rV[gg + nrW * aa];
	}
      }
    } else {  /* ncA == 1 */ /* nc B == 1 */
      for (i=0; i<nrA; i++){
	aa = (int) rA[i]-1;
	for (j=0; j<nrB; j++){
	  //gg = (int) rB[j + nrB]-1;
	  gg = (int) rB[j]-1;
	  //Rprintf("i %i j %i aa: %i gg: %i rW-idx %i rV-idx %i \n",
	  //  i,j, aa, gg, aa+nrV*gg, gg+nrW*aa);
	  ans = ans +	
	    rW[(aa + nrV * gg)] * rV[(gg + nrW * aa)];
	}
      }
    }
  }

  return(ans);
}



//[[Rcpp::export]]
NumericVector trAWBlist_(List Alist, NumericMatrix W, List Blist, int mode=0)
{

  int nA = Alist.length(), nB = Blist.length();
  int Astart = 0, idx=0;

  NumericVector ans(nA * nB);

  for (int Bjj=0; Bjj < nB; Bjj++){
    NumericMatrix Bitem = Blist[Bjj];
    // Rf_PrintValue(Bitem);
    if (mode != 0)
      Astart = Bjj;
    
    for (int Aii=Astart; Aii < nA; Aii++){
      NumericMatrix Aitem = Alist[Aii];
      // Rf_PrintValue(Aitem);
      ans[idx++] = trAWB_(Aitem, W, Bitem);
    }
  }
  
  return(ans);
}

//[[Rcpp::export]]
NumericVector trAWBWlist_(List Alist, NumericMatrix W, List Blist, int mode=0)
{

  int nA = Alist.length(), nB = Blist.length();
  int Astart = 0, idx=0;

  NumericVector ans(nA * nB);

  for (int Bjj=0; Bjj < nB; Bjj++){
    NumericMatrix Bitem = Blist[Bjj];
    // Rf_PrintValue(Bitem);
    if (mode != 0)
      Astart = Bjj;
    
    for (int Aii=Astart; Aii < nA; Aii++){
      NumericMatrix Aitem = Alist[Aii];
      // Rf_PrintValue(Aitem);
      ans[idx++] = trAWBW_(Aitem, W, Bitem);
    }
  }
  
  return(ans);
}




























// //[[Rcpp::export]]
// void foo(NumericMatrix A){
//   Rcout << "hej" << std::endl;
//   for (int i=0; i<4; i++){
//     double e = A[i];
//     Rcout << "i " << i << "e " << e << std::endl;
//   }
// }


// double trAWprim_(NumericMatrix A, NumericMatrix W)
// {

//   int nrA = A.nrow(), ncA = A.ncol();
//   int nrW = W.nrow(), ncW = W.ncol();

//   double ans=0;
//   int ii, idx;

//   if (ncA == 2){
//     for (ii=0; ii<nrA; ii++){
//       idx = (int) (A[ii] - 1 + nrW * (A[ii + nrA] - 1));
//       //Rprintf("%i ", idx);
//       ans = ans + W[idx];
//     }
//     ans = 2 * ans;    
//   } else {
//     for (ii=0; ii < nrA; ii++){
//       idx = (int) ((A[ii] - 1) * ( nrW + 1));
//       //Rprintf("%i %f \n", idx, rA[ii]);
//       ans = ans + W[idx];
//     }
//   }
//   return(ans);
// }

















// //[[Rcpp::export]]
// double trAW_(NumericMatrix A, NumericMatrix W)
// {
//   //double *rans;
//   // SEXP Adims, Wdims, ans;

//   int nrA = A.nrow(), ncA = A.ncol();
//   int nrW = W.nrow(), ncW = W.ncol();

//   // IntegerVector rA = (IntegerVector) A;
//   // IntegerVector rW = (IntegerVector) W;

//   double *rA, *rW;
//   rA = REAL(A);
//   rW = REAL(W);

//   SEXP ans;
//   PROTECT(ans =Rf_allocVector(REALSXP,1));
//   double *rans;
//   rans      = REAL(ans);
//   *rans = 1000.0;
//   trAWprim_(rA, &nrA, &ncA, rW, &nrW, &ncW, ans);
  
//   // UNPROTECT(1);
//   return(ans);
// }



  // PROTECT(ans =allocVector(REALSXP,1));
  // rans      = REAL(ans);
  // *rans = 1000.0;



  // trAWprim(rA, &nrA, &ncA,
  // 	   rW, &nrW, &ncW, rans);

  //Rprintf("trAWB: %f\n",rans);

  
  // Adims = getAttrib(A, R_DimSymbol);
  // Wdims = getAttrib(W, R_DimSymbol);

  // PROTECT(A = AS_NUMERIC(A));
  // PROTECT(W = AS_NUMERIC(W));


  // nrA = INTEGER(Adims)[0];  ncA = INTEGER(Adims)[1];
  // nrW = INTEGER(Wdims)[0];  ncW = INTEGER(Wdims)[1];
