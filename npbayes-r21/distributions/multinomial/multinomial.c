#ifndef MULTINOMIAL
#define MULTINOMIAL

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include "../../utilities/mxutils.c"

typedef struct {
  int numdim;
  double *eta;
} HH;

typedef int SS;

typedef int* QQ;

int numitems(QQ qq) {
  return qq[0];
}

double marglikelihood(HH hh, QQ qq, SS ss) {
  return log( (hh.eta[ss]+qq[ss])/(hh.eta[0]+qq[0]) );
}

double marglikelihoods(double *clik, HH hh, int numqq, QQ *qq, SS ss) {
  double etas, eta0;
  int ii;
  eta0 = hh.eta[0];
  etas = hh.eta[ss];
  for ( ii = 0 ; ii < numqq ; ii++ ) 
    clik[ii] = (etas + qq[ii][ss]) / (eta0 + qq[ii][0]);
  return 0;
}

void adddata(HH hh, QQ qq, SS ss) {
  qq[ss] ++;
  qq[0] ++;
}

double adddatalik(HH hh, QQ qq, SS ss) {
  return log( (hh.eta[ss] + (qq[ss] ++)) / (hh.eta[0] + (qq[0] ++)) );
}

void deldata(HH hh, QQ qq, SS ss) {
  qq[ss] --;
  qq[0] --;
}

QQ newclass(HH hh) {
  int ii;
  QQ result;
  result = mxMalloc(sizeof(int)*(hh.numdim+1));
  for ( ii = 0 ; ii <= hh.numdim ; ii++ ) 
    result[ii] = 0;
  return result;
}

SS *mxReadSSVector(const mxArray *mvector) {
  SS *result;
  double *pr;
  int ii;
  pr = mxGetPr(mvector);
  result = mxMalloc(sizeof(SS)*mxGetNumberOfElements(mvector));
  for ( ii = 0 ; ii < mxGetNumberOfElements(mvector) ; ii++ )
    result[ii] = pr[ii];
  return result;
}

mxArray *mxWriteSSVector(int numss, SS *ss) {
  mxArray *result;
  double *pr;
  int ii;
  result = mxCreateDoubleMatrix(1,numss, mxREAL);
  pr = mxGetPr(result);
  for ( ii = 0 ; ii < numss ; ii++ )
    pr[ii] = ss[ii];
  mxFree(ss);
  return result;
}

void mxFreeSSVector(int numss, SS *ss) {
  mxFree(ss);
}

HH mxReadHH(const mxArray *mvector) {
  double *pr, sum;
  int ii, jj;
  HH result;
  result.numdim = mxGetNumberOfElements(mvector);
  result.eta    = mxMalloc(sizeof(double)*(1+result.numdim));
  sum           = 0.0;
  pr            = mxGetPr(mvector);
  for ( ii = 0 ; ii < result.numdim ; ii++ ) 
    sum += result.eta[ii+1]  = pr[ii];
  result.eta[0] = sum;
  return result;
}

void mxFreeHH(HH hh) {
  mxFree(hh.eta);
}

QQ *mxReadQQVector(HH hh, const mxArray *marray, int maxnum) {
  double *pr;
  int ii, jj, mm, nn, sum;
  QQ *result;
  mm = mxGetM(marray);
  nn = mxGetN(marray);
  if ( mm != hh.numdim ) mexErrMsgTxt("Number of dimensions don't match.");
  maxnum = max(maxnum,nn);
  result = mxMalloc(sizeof(QQ)*maxnum);
  pr = mxGetPr(marray);
  for ( jj = 0 ; jj < nn ; jj++ ) {
    result[jj] = newclass(hh);
    sum = 0;
    for ( ii = 0 ; ii < mm ; ii++ ) {
      sum += result[jj][ii+1] = pr[ii+jj*mm];
    }
    result[jj][0] = sum;
  }
  for ( jj = nn ; jj < maxnum ; jj++ ) 
    result[jj] = newclass(hh);
  return result;
}

mxArray *mxWriteQQVector(HH hh, int nn, int maxnum, QQ *qq) {
  mxArray *result;
  double *pr;
  int ii,jj;
  result = mxCreateDoubleMatrix(hh.numdim,nn,mxREAL);
  pr = mxGetPr(result);
  mxdebug2(3,"Write qq (%dx%d): ",hh.numdim,nn);
  for ( jj = 0 ; jj < nn ; jj ++ ) {
    for ( ii = 0 ; ii < hh.numdim ; ii++ )
      pr[ii+jj*hh.numdim] = qq[jj][ii+1];
    mxFree(qq[jj]);
    mxdebug1(4,"%d ",jj);
  }
  mxdebug0(3,"Free qq: ");
  for ( jj = nn ; jj < maxnum ; jj ++ ) {
    mxFree(qq[jj]);
    mxdebug1(4,"%d ",jj);
  }
  mxFree(qq);
  mxdebug0(3,"\n");
  return result;
}

#endif
