#include "../distributions/multinomial/multinomial.c"
#include "hdp.c"

#define HDPIN prhs[0]
#define NUMITER prhs[1]
#define DOCONPARAM prhs[2]
#define DOLIK prhs[3]
#define DODEBUG prhs[4]
#define HDPOUT plhs[0]
#define LIK plhs[1]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
  int ii, jj, cc, ss;
  mxArray *tmparray, *tmpcell;
  double *tmpdouble;
  HDP *hdp;

  if (nrhs!=4 && nrhs!=5) mexErrMsgTxt("Four/five input argument required.");
  if (!mxIsStruct(HDPIN)) mexErrMsgTxt("HDP structure expected.");
  if (mxGetNumberOfElements(HDPIN)!=1) mexErrMsgTxt("One structure expected.");

#ifndef NODEBUG
  DEBUG = nrhs==4 ? 0 : *mxGetPr(DODEBUG);
#endif

  hdp = mxReadHDP(HDPIN);

  LIK = mxCreateDoubleMatrix(1, *mxGetPr(NUMITER),mxREAL);

  mxdebug0(1,"Running hdpMultinomial_iterate.\n");
  hdp_iterate(hdp, mxGetPr(LIK), 
      *mxGetPr(NUMITER), *mxGetPr(DOCONPARAM), *mxGetPr(DOLIK));
  mxdebug0(1,"Done hdpMultinomial_iterate.\n");

  HDPOUT = mxDuplicateArray(HDPIN);
  mxWriteHDP(HDPOUT,hdp);
}


