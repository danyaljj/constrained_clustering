#ifndef MXUTILS
#define MXUTILS
#include "mex.h"
#include <math.h>

#define max(x1,x2) ( (x1) < (x2) ? (x2) : (x1) )

#ifndef NODEBUG
int DEBUG;
#define mxdebug0(num,string)       if(DEBUG>=(num))mexPrintf(string);
#define mxdebug1(num,string,a)     if(DEBUG>=(num))mexPrintf(string,a);
#define mxdebug2(num,string,a,b)   if(DEBUG>=(num))mexPrintf(string,a,b);
#define mxdebug3(num,string,a,b,c) if(DEBUG>=(num))mexPrintf(string,a,b,c);
#define mxdebug4(num,string,a,b,c,d) if(DEBUG>=(num))mexPrintf(string,a,b,c,d);
#define mxdebugarray(num,string,str,array,length) { \
  if (DEBUG >= (num)) { \
    int ii; \
    mexPrintf("%s: ",string); \
    for ( ii = 0 ; ii < length ; ii++) { \
      mexPrintf(str,array[ii]); \
      if (ii % 5 == 4) { \
        mexPrintf(" : "); \
      } else { \
        mexPrintf(" ");  \
      } \
    } \
    mexPrintf("\n"); \
  } \
}
#else
#define mxdebug0(num,string)     
#define mxdebug1(num,string,a)     
#define mxdebug2(num,string,a,b)  
#define mxdebug3(num,string,a,b,c)
#define mxdebug4(num,string,a,b,c,d)
#define mxdebugarray(num,string,str,array,length) 
#endif

mxArray *mxReadField(const mxArray *mstruct,const char *fieldstr) {
  mxArray *result;
  result = mxGetField(mstruct,0,fieldstr); 
  if (result == NULL) { 
    mexPrintf(fieldstr); 
    mexErrMsgTxt(" field missing."); 
  } 
  if (DEBUG>=3) mexPrintf("Read %s.\n",fieldstr); 
  return result;
}

void mxWriteField(mxArray *mstruct,const char *fieldstr,mxArray *mfield) {
  if ( mxGetField(mstruct,0,fieldstr) != NULL ) 
    mxDestroyArray(mxGetField(mstruct,0,fieldstr)); 
  mxSetField(mstruct,0,fieldstr,mfield);
  if (DEBUG>=3) mexPrintf("Write %s.\n",fieldstr); 
}


double mxReadScalar(const mxArray *mscalar) { 
  if (DEBUG>=3) mexPrintf("Value = %g.\n",*mxGetPr(mscalar));
  return (*mxGetPr(mscalar));
}

mxArray *mxWriteScalar(double var) { 
  mxArray *result;
  if (DEBUG>=3) mexPrintf("Value = %g.\n",var);
  result = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(result) = var; 
  return result;
}

#define mxReadVectorDef(funcname,type,str) \
  type *funcname(const mxArray *mvector,int number,type shift,type init) { \
    double *mdouble; \
    type *result;  \
    int ii; \
    number = max(number,mxGetNumberOfElements(mvector)); \
    result = mxMalloc(sizeof(type)*number); \
    mdouble = mxGetPr(mvector); \
    for ( ii = 0 ; ii < mxGetNumberOfElements(mvector) ; ii++ ) \
      result[ii] = mdouble[ii] + shift; \
    for ( ii = mxGetNumberOfElements(mvector) ; ii < number ; ii++ ) \
      result[ii] = init; \
    if (DEBUG>=3) { \
      mexPrintf("Value = "); \
      for ( ii = 0 ; ii < number ; ii++ ) \
        mexPrintf(str,result[ii]); \
      mexPrintf("\n"); \
    } \
    return result; \
  } 
mxReadVectorDef(mxReadIntVector,int,"%d ");
mxReadVectorDef(mxReadDoubleVector,double,"%g ");

#define mxWriteVectorDef(funcname, type, str) \
  mxArray *funcname(int mm,int nn,type *var,type shift) { \
    mxArray *result; \
    double *mdouble; \
    int ii; \
    result = mxCreateDoubleMatrix(mm,nn,mxREAL); \
    mdouble = mxGetPr(result); \
    for ( ii = 0 ; ii < mm*nn ; ii++ ) \
      mdouble[ii] = var[ii] + shift; \
    if (DEBUG>=3) { \
      mexPrintf("Value = "); \
      for ( ii = 0 ; ii < mm*nn ; ii++ ) \
        mexPrintf(str,var[ii]); \
      mexPrintf("\n"); \
    } \
    mxFree(var); \
    return result; \
  } 
mxWriteVectorDef(mxWriteIntVector, int, "%d ");
mxWriteVectorDef(mxWriteDoubleVector, double, "%g ");


#define mxReadCellVectorDef(funcname,type) \
  type **funcname(const mxArray *mcell, type shift) { \
    mxArray *mvector; \
    double *mdouble; \
    int ii, jj; \
    type **result; \
    result = mxMalloc(sizeof(type*)*mxGetNumberOfElements(mcell));  \
    for ( jj = 0 ; jj < mxGetNumberOfElements(mcell) ; jj++ ) {  \
      mvector = mxGetCell(mcell,jj);  \
      mdouble = mxGetPr(mvector);  \
      result[jj] = mxMalloc(sizeof(type)*mxGetNumberOfElements(mvector));  \
      for ( ii = 0 ; ii < mxGetNumberOfElements(mvector) ; ii++ )  \
        result[jj][ii] = mdouble[ii] + shift;  \
    } \
    return result; \
  }
mxReadCellVectorDef(mxReadIntCellVector,int);
mxReadCellVectorDef(mxReadDoubleCellVector,double);

#define mxWriteCellVectorDef(funcname, type) \
  mxArray *funcname(int numcell,int *numentry,type **var,type shift) { \
    mxArray *result, *mvector; \
    double *mdouble; \
    int ii, jj; \
    result = mxCreateCellMatrix(1,numcell); \
    for ( jj = 0 ; jj < numcell ; jj++) { \
      mvector = mxCreateDoubleMatrix(1,numentry[jj],mxREAL); \
      mxSetCell(result,jj,mvector); \
      mdouble = mxGetPr(mvector); \
      for ( ii = 0 ; ii < numentry[jj] ; ii++ ) \
        mdouble[ii] = var[jj][ii] + shift; \
      mxFree(var[jj]); \
    } \
    mxFree(var); \
    return result; \
  }
mxWriteCellVectorDef(mxWriteIntCellVector, int);
mxWriteCellVectorDef(mxWriteDoubleCellVector, double);

#define mxReadMatrixDef(funcname,type) \
  type **funcname(const mxArray *marray,int mm,int nn,type shift,type init) { \
    double *mdouble; \
    int ii, jj, m1, n1; \
    type **result; \
    mdouble = mxGetPr(marray); \
    mm = max(mm, m1 = mxGetM(marray)); \
    nn = max(nn, n1 = mxGetN(marray)); \
    result = mxMalloc(sizeof(type*)*mm); \
    for ( jj = 0 ; jj < mm ; jj++ ) { \
      result[jj] = mxMalloc(sizeof(type)*nn); \
    } \
    for ( jj = 0 ; jj < m1 ; jj++ ) {\
      for ( ii = 0 ; ii < n1 ; ii++ ) \
        result[jj][ii] = mdouble[ii*m1+jj] + shift; \
      for ( ii = n1 ; ii < nn ; ii++ ) \
        result[jj][ii] = init; \
    } \
    for ( jj = m1 ; jj < mm ; jj++ ) \
      for ( ii = 0 ; ii < nn ; ii++ ) \
        result[jj][ii] = init; \
    return result; \
  }
mxReadMatrixDef(mxReadIntMatrix,int);
mxReadMatrixDef(mxReadDoubleMatrix,double);

#define mxWriteMatrixDef(funcname, type) \
  mxArray *funcname(int mm,int nn,int maxm,type **var,type shift) { \
    mxArray *result; \
    double *mdouble; \
    int ii, jj; \
    result  = mxCreateDoubleMatrix(mm,nn,mxREAL); \
    mdouble = mxGetPr(result); \
    for ( jj = 0 ; jj < mm ; jj++) { \
      for ( ii = 0 ; ii < nn ; ii++ ) \
        mdouble[jj+mm*ii] = var[jj][ii] + shift; \
      mxFree(var[jj]); \
    } \
    for ( jj = mm ; jj < maxm ; jj ++ ) \
      mxFree(var[jj]); \
    mxFree(var); \
  }
mxWriteMatrixDef(mxWriteIntMatrix, int);
mxWriteMatrixDef(mxWriteDoubleMatrix, double);


void mxFreeCellVector(void **var,int numcell) {
  int jj; 
  for ( jj = 0 ; jj < numcell ; jj++ )  
    mxFree(var[jj]); 
  mxFree(var); 
}

#endif
