#ifndef _BASE
#define _BASE
/*
 * Structure representing the base distribution, in particular, the 
 * parameters for the base distribution, and samples from the base 
 * distribution (one sample per class).
 *
 * numclass     The number of samples from the base distribution (number of
 *              clusters in HDP mixture).
 * maxclass     The number of samples for which we have allocated memory.
 *              This is initialized to (numclass+2)*2.
 * hh           The parameters for the base distribution.
 * classqq      The sufficient statistics for data items associated with
 *              each cluster.  In matlab, classqq(:,kk) refers to the
 *              statistics for class kk.  There is numclass+1 represented
 *              clusters (one additional cluster with no data associated 
 *              with it).
 * beta         Only used internally.  A vector of numclass 0's following by
 *              one 1.  In matlab this is a row vector.
 *
 * BASE *mxReadBase(mxArray *mstruct);
 *              Reads in a BASE struct from a matlab struct.
 * void mxWriteBase(mxArray *result,BASE *base);
 *              Writes a BASE struct to a matlab struct.  Overwrites fields
 *              in result if necessary.  Frees memory allocated.
 */

#include "../utilities/mxutils.c"

typedef struct {
  int numclass, maxclass;
  HH hh;
  QQ *classqq;
  double *beta;
} BASE;

BASE *mxReadBase(mxArray *mstruct) {
  BASE *result;
  int ii, maxclass;
  result = mxMalloc(sizeof(BASE));
  result->numclass = mxReadScalar(mxReadField(mstruct,"numclass"));
  result->maxclass = maxclass = (result->numclass+2) * 2;
  result->hh       = mxReadHH(mxReadField(mstruct,"hh"));
  result->classqq  = mxReadQQVector(result->hh,mxReadField(mstruct,"classqq"),
                        result->maxclass);
  result->beta = mxMalloc(sizeof(double)*maxclass);
  for ( ii = 0 ; ii < maxclass ; ii++) 
    result->beta[ii] = 0.0;
  result->beta[result->numclass] = 1.0;
  return result;
}

void mxWriteBase(mxArray *result,BASE *base) {
  mxWriteField(result,"numclass",mxWriteScalar(base->numclass));
  mxWriteField(result,"classqq",mxWriteQQVector(base->hh,
        base->numclass+1,base->maxclass,base->classqq));
  mxFreeHH(base->hh);
  mxFree(base->beta);
  mxFree(base);
}

#endif
