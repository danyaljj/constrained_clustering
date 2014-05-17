#include "../distributions/multinomial/multinomial.c"
#include "hdp.c"

#define HDPIN prhs[0]
#define NUMBURNIN prhs[1]
#define NUMSAMPLE prhs[2]
#define DOCONPARAM prhs[3]
#define PREDICTJJ prhs[4]
#define DODEBUG prhs[5]
#define HDPOUT plhs[0]
#define LIK plhs[1]

void mexFunction( int nlhs, mxArray *plhs[], 
                          int nrhs, const mxArray*prhs[] )
     
{
  int numpredict, *predictjj;
  int ii, jj, cc, ss;
  mxArray *tmparray, *tmpcell;
  double *tmpdouble;
  HDP *hdp;

  if (nrhs!=5 && nrhs!=6) mexErrMsgTxt("Four/five input argument required.");
  if (!mxIsStruct(HDPIN)) mexErrMsgTxt("HDP structure expected.");
  if (mxGetNumberOfElements(HDPIN)!=1) mexErrMsgTxt("One structure expected.");
#ifndef NODEBUG
  DEBUG = nrhs==5 ? 0 : *mxGetPr(DODEBUG);
#endif

  hdp = mxReadHDP(HDPIN);

  numpredict = mxGetNumberOfElements(PREDICTJJ);
  tmpdouble  = mxGetPr(PREDICTJJ);
  predictjj  = mxMalloc(sizeof(int)*numpredict);
  for ( jj = 0 ; jj < numpredict ; jj++ )
    predictjj[jj] = ((int)tmpdouble[jj]) - 1;

  LIK = mxCreateDoubleMatrix(*mxGetPr(NUMSAMPLE),numpredict,mxREAL);

  mxdebug0(1,"Running hdp_predict.\n");
  hdp_predict(hdp, mxGetPr(LIK), *mxGetPr(NUMBURNIN), *mxGetPr(NUMSAMPLE),
      numpredict, predictjj, *mxGetPr(DOCONPARAM));
  mxdebug0(1,"Done hdp_predict.\n");

  HDPOUT = mxDuplicateArray(HDPIN);
  mxWriteHDP(HDPOUT,hdp);
}


