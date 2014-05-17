#ifndef _CONPARAM
#define _CONPARAM
/*
 * Structure representing a concentration parameter in the HDP.  In
 * particular, we represent the parameter itself, and parameters for
 * a gamma prior over it, the number of DPs using the concentration
 * parameter, and the total number of data items and tables for each DP.
 *
 * alpha        The concentration parameter itself.  Has to be nonnegative.
 * alphaa       The shape parameter for the gamma prior over alpha.
 * alphab       The inverse scale parameter for the gamma prior over alpha.
 * numdp        The number of DPs using this concentration parameter.
 * totalnd      The number of data items in each DP.  In matlab this is a row
 *              vector with jj'th entry belonging to the jj'th DP associated
 *              this concentration parameter.
 * totalnt      Same as totalnd for number of tables in each DP.
 *
 * CONPARAM *mxReadConparamVector(mxArray *mcell);
 *              reads in a CONPARAM struct from a matlab struct. 
 * void mxWriteConparamVector(mxArray *result,CONPARAM *cparray);
 *              writes a CONPARAM struct to a matlab struct.  Overwrites 
 *              fields if necessary.  Frees memory allocated.
 */

#include "../utilities/mxutils.c"

typedef struct {
  double alpha, alphaa, alphab;
  int numdp, *totalnd, *totalnt;
} CONPARAM;

CONPARAM *mxReadConparamVector(mxArray *mcell) {
  mxArray *mstruct;
  CONPARAM *result, *cp;
  int ii, nn;
  nn = mxGetNumberOfElements(mcell);
  result = mxMalloc(sizeof(CONPARAM)*nn);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    mstruct = mxGetCell(mcell,ii);
    cp  = result + ii;
    cp->alpha = mxReadScalar(mxReadField(mstruct,"alpha"));
    cp->alphaa = mxReadScalar(mxReadField(mstruct,"alphaa"));
    cp->alphab = mxReadScalar(mxReadField(mstruct,"alphab"));
    cp->numdp  = mxReadScalar(mxReadField(mstruct,"numdp"));
    cp->totalnd = mxReadIntVector(mxReadField(mstruct,"totalnd"),0,0,0);
    cp->totalnt = mxReadIntVector(mxReadField(mstruct,"totalnt"),0,0,0);
  }
  return result;
}

void mxWriteConparamVector(mxArray *result,CONPARAM *cparray) {
  mxArray *mstruct;
  CONPARAM *cp;
  int ii, nn;
  nn = mxGetNumberOfElements(result);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    cp = cparray + ii;
    mstruct = mxGetCell(result,ii);
    mxWriteField(mstruct,"alpha",mxWriteScalar(cp->alpha));
    mxWriteField(mstruct,"totalnd",mxWriteIntVector(1,cp->numdp,cp->totalnd,0));
    mxWriteField(mstruct,"totalnt",mxWriteIntVector(1,cp->numdp,cp->totalnt,0));
  }
  mxFree(cparray);
}

#endif
