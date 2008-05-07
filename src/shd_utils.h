
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
