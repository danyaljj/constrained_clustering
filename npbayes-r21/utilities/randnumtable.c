#include "randutils.c"

#define mxIN1 prhs[0]
#define mxIN2 prhs[1]
#define mxOUT plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
                      int nrhs, const mxArray*prhs[] )
     
{
  double *xx, *xxend, *al, *nd;

  if (nrhs!=2) mexErrMsgTxt("Two input argument required.");
  if (!mxIsDouble(mxIN1) || !mxIsDouble(mxIN2)) 
    mexErrMsgTxt("Double array expected.");
  if (mxGetNumberOfElements(mxIN1) != mxGetNumberOfElements(mxIN2))
    mexErrMsgTxt("Array sizes must be equal.");


  mxOUT = mxDuplicateArray(mxIN2);
  xx = mxGetPr(mxOUT);
  al = mxGetPr(mxIN1);
  nd = mxGetPr(mxIN2);
  xxend = xx + mxGetNumberOfElements(mxOUT);
  while ( xx < xxend ) {
    *xx++ = (double) randnumtable(*al++,*nd++);
  }
}



