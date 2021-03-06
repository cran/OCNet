/* Reinhard Furrer,  fall 2019, based on spam's version and 'Writing R extensions'. */

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>

#include "invperm.c"

/* to get all functions:
   nm -g ./OCNet.so | grep " T "
*/



/* .Call calls */

static const R_CallMethodDef CallEntries[] = {
  {"inv_permutation", (DL_FUNC) &inv_permutation, 1},
  {NULL, NULL, 0}
};


/* .Fortran calls */

extern void F77_NAME(allinone    )( void *, void *, void *, void *, void *, void *,void *, void *,  void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(permandsolve)( void *, void *, void *, void *, void *, void *,void *, void *,  void *, void *, void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
    {"permandsolve",        (DL_FUNC) &F77_NAME(permandsolve ),13},
    {"allinone",            (DL_FUNC) &F77_NAME(allinone     ),21},
    {NULL, NULL, 0}
};

void R_init_OCNet(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

