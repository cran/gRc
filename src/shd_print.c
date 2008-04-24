#include <Rdefines.h>
#include <R.h>

/* Forward declarations */

void printveci(int *x, int *n);
void printvecd(double *x, int *n);
void printmati(int *Avec, int *nr, int *nc);
void printmatd(double *Avec, int *nr, int *nc);


void printveci(int *x, int *n){

  for (int i=0; i<*n; i++)
    Rprintf(" %2i ",x[i]);
  Rprintf("\n");
}

void printvecd(double *x, int *n){

  for (int i=0; i<*n; i++)
    Rprintf(" %f ",x[i]);
  Rprintf("\n");
}

void printmati(int *Avec, int *nr, int *nc){
  int ii, jj;

  for (ii=0; ii<*nr; ii++){
    for (jj=0; jj<*nc; jj++){
      Rprintf(" %i ", (int) Avec[ii + *nr * jj]);
    }
    Rprintf("\n");
  }
  Rprintf(" ---------------------------\n");
}

void printmatd(double *Avec, int *nr, int *nc){
  int ii, jj;
  
  for (ii=0; ii<*nr; ii++){
    for (jj=0; jj<*nc; jj++){
      Rprintf(" %10.5f ", Avec[ii + *nr * jj]);
    }
    Rprintf("\n");
  }
  Rprintf(" ---------------------------\n");
}
