#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _gRc_RcppExport_registerCCallable(void);
extern SEXP _gRc_trAWBV_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRc_trAWBW_(SEXP, SEXP, SEXP);
extern SEXP _gRc_trAWBWlist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRc_trAWB_(SEXP, SEXP, SEXP);
extern SEXP _gRc_trAWBlist_(SEXP, SEXP, SEXP, SEXP);
extern SEXP _gRc_trAW_(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_gRc_RcppExport_registerCCallable", (DL_FUNC) &_gRc_RcppExport_registerCCallable, 0},
    {"_gRc_trAWBV_",                      (DL_FUNC) &_gRc_trAWBV_,                      4},
    {"_gRc_trAWBW_",                      (DL_FUNC) &_gRc_trAWBW_,                      3},
    {"_gRc_trAWBWlist_",                  (DL_FUNC) &_gRc_trAWBWlist_,                  4},
    {"_gRc_trAWB_",                       (DL_FUNC) &_gRc_trAWB_,                       3},
    {"_gRc_trAWBlist_",                   (DL_FUNC) &_gRc_trAWBlist_,                   4},
    {"_gRc_trAW_",                        (DL_FUNC) &_gRc_trAW_,                        2},
    {NULL, NULL, 0}
};

void R_init_gRc(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
