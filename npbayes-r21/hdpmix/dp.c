#ifndef _DP
#define _DP

/*
 * Structure representing each Dirichlet process.  In particular, the
 * concentration parameter, the parent distribution, the beta weights, 
 * the numbers of data items/tables, and the statistics of data items 
 * directly associated with this DP.
 *
 * alpha        The concentration parameter itself.
 * numdata      The number of data items drawn directly from this DP.  This
 *              excludes virtual data items drawn as a result of data items
 *              associated with children DPs.
 * datacc       The cluster to which data items are associated with.  In
 *              matlab this is a row vector and datacc(ii) is the cluster of
 *              data item ii.
 * datass       The sufficient statistics of data items.  In matlab
 *              datass(:,ii) is the statistics for data item ii.
 * classnd      The number of data items in each component, including virtual
 *              data items.  There are numclass+1 entries, one for each
 *              component and an extra one for the empty class.  In matlab 
 *              this is a row vector, with cc'th entry being the number of 
 *              data items in component cc.
 * classnt      Same as classnd, except for number of tables.
 * beta         The weight associated with each component in this DP.  There
 *              are numclass+1 entries, one for each component, plus an
 *              additional one for the empty class.
 *
 * DP *mxReadDPVector(mxArray *mcell,int *dpstate,int mm)
 *              Reads in a DP struct from a matlab struct.  Allocates enough
 *              memory for mm components.
 * void mxWriteDPVector(mxArray *result,int nn,DP *dparray,int *dpstate,
 *                      int numclass)
 *              Writes a DP struct to a matlab struct.  Overwrites fields in
 *              result if necessary.  Frees memory allocated.
 */

#include "../utilities/mxutils.c"

#define ACTIVE 2
#define FROZEN 1
#define HELDOUT 0

typedef struct DP{
  double alpha, *beta;
  int *classnd, *classnt;
  int numdata, *datacc;
  SS *datass;
} DP;

DP *mxReadDPVector(mxArray *mcell, int *dpstate, int mm) {
  mxArray *mstruct;
  DP *result, *dp;
  int ii, nn, pp;
  nn = mxGetNumberOfElements(mcell);
  result = mxMalloc(sizeof(DP)*nn);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    mstruct       = mxGetCell(mcell,ii);
    dp            = result+ii;
    dp->numdata   = mxReadScalar(mxReadField(mstruct,"numdata"));
    dp->datass    = mxReadSSVector(mxReadField(mstruct,"datass"));
    if ( dpstate[ii] == ACTIVE || dpstate[ii] == FROZEN ) {
      int cc;
      dp->classnd = mxReadIntVector(mxReadField(mstruct,"classnd"),mm,0,0);
      dp->classnt = mxReadIntVector(mxReadField(mstruct,"classnt"),mm,0,0);
      dp->datacc  = mxReadIntVector(mxReadField(mstruct,"datacc"),0,-1,0);
    } else {
      dp->classnd = NULL;
      dp->classnt = NULL;
      dp->datacc  = NULL;
    }
    if ( dpstate[ii] == ACTIVE ) {
      dp->alpha   = mxReadScalar(mxReadField(mstruct,"alpha"));
      dp->beta    = mxReadDoubleVector(mxReadField(mstruct,"beta"),mm,0.0,0.0);
    } else {
      dp->alpha   = 0.0;
      dp->beta    = NULL;
    }
  }
  return result;
}

void mxWriteDPVector(mxArray *result,int nn,DP *dparray,int *dpstate,
    int numclass) {
  mxArray *mstruct;
  DP *dp;
  int ii;
  mxdebug2(3,"Write DP %d-vector, numclass=%d\n",nn,numclass);
  for ( ii = 0 ; ii < nn ; ii++ ) {
    dp = dparray+ii;
    mstruct = mxGetCell(result,ii);
    if ( dpstate[ii] == ACTIVE ) {
      mxWriteField(mstruct,"alpha",mxWriteScalar(dp->alpha));
      mxWriteField(mstruct,"beta",mxWriteDoubleVector(1,numclass,dp->beta,0.0));
      mxFreeSSVector(dp->numdata,dp->datass);
    } else {
      mxWriteField(mstruct,"alpha",mxCreateDoubleMatrix(0,0,mxREAL));
      mxWriteField(mstruct,"beta",mxCreateDoubleMatrix(0,0,mxREAL));
    }

    if ( dpstate[ii] == ACTIVE || dpstate[ii] == FROZEN ) {
      int cc;
      mxWriteField(mstruct,"classnd",
          mxWriteIntVector(1,numclass,dp->classnd,0));
      mxWriteField(mstruct,"classnt",
          mxWriteIntVector(1,numclass,dp->classnt,0));
      mxWriteField(mstruct,"datacc",
          mxWriteIntVector(1,dp->numdata,dp->datacc,1));
    } else {
      mxWriteField(mstruct,"classnd",mxCreateDoubleMatrix(0,0,mxREAL));
      mxWriteField(mstruct,"classnt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxWriteField(mstruct,"datacc",mxCreateDoubleMatrix(0,0,mxREAL));
    }
  }
  mxFree(dparray);
}

#endif
