#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP rconipm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP scorematching(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP trAW(SEXP, SEXP);
extern SEXP trAWB(SEXP, SEXP, SEXP);
extern SEXP trAWBlist(SEXP, SEXP, SEXP, SEXP);
extern SEXP trAWBV(SEXP, SEXP, SEXP, SEXP);
extern SEXP trAWBW(SEXP, SEXP, SEXP);
extern SEXP trAWBWlist(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rconipm",       (DL_FUNC) &rconipm,       10},
    {"scorematching", (DL_FUNC) &scorematching,  6},
    {"trAW",          (DL_FUNC) &trAW,           2},
    {"trAWB",         (DL_FUNC) &trAWB,          3},
    {"trAWBlist",     (DL_FUNC) &trAWBlist,      4},
    {"trAWBV",        (DL_FUNC) &trAWBV,         4},
    {"trAWBW",        (DL_FUNC) &trAWBW,         3},
    {"trAWBWlist",    (DL_FUNC) &trAWBWlist,     4},
    {NULL, NULL, 0}
};

void R_init_gRc(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
