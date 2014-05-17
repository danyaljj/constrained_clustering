#include "randutils.c"

#define mxIN prhs[0]
#define mxOUT plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
                      int nrhs, const mxArray*prhs[] )
     
{
  double *xx, *xxend;

  if (nrhs!=1) mexErrMsgTxt("One input argument required.");
  if (!mxIsDouble(mxIN)) mexErrMsgTxt("Double array expected.");

  mxOUT = mxDuplicateArray(mxIN);
  xx = mxGetPr(mxOUT);
  for (xxend = xx + mxGetNumberOfElements(mxOUT) ; xx < xxend ; xx++) {
    *xx = randgamma(*xx);
  }
}



