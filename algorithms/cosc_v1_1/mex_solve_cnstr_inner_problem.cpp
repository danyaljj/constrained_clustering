//
// (C)2012 Syama Sundar Rangapuram and Matthias Hein
// Max Planck Institute for Computer Science, Saarbruecken
// Machine Learning Group, Saarland University, Germany
// http://www.ml.uni-saarland.de
//    

#include <math.h>
#include "mex.h"
#include "matrix.h"

// double* proj_smplx( double*, double );
double* ProjKSimplex(double* f,int len, int K);
double BacktrackLineSearch( double*, double *alpha, double *v, double *D, double* Xobs, double* r2, mwIndex *irs, mwIndex *jcs, double L, int n, int m, int cols, double C, double, double, double, double* );

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[]) {
             
//   printf("start of the C code\n");
//   printf("%d, %d\n", nrhs, nlhs);
  // Test number of parameters.
  if (nrhs != 11 || nlhs != 5) {
    mexWarnMsgTxt("Usage: [X,rval,v,Obj,niter]=solveInnerProblem(W,s,alpha,L,r2,u,C, MAXITER,EPS)");
    return;
  }
  // get important parameters
  int rows = (int)mxGetM(prhs[0]); // number of rows of W
  int cols = (int)mxGetN(prhs[0]); // number of columns of W (should be the same)
  int len  = (int)mxGetM(prhs[1]); // the desired output
  int lenrval = (int)mxGetM(prhs[2]); // rval

  if(!mxIsSparse(prhs[0])) { mexWarnMsgTxt("Matrix is not sparse");}
 
  if(rows!=cols){
    mexWarnMsgTxt("Sparse matrix is not square");
    return;
  }

  if(rows!=len){
    mexWarnMsgTxt("Length of the vector is not the same as the number of the rows of the sparse matrix");
    return;
  }
  
  //printf("the correct function");
    
  // Create output array and compute values
  double* sr = mxGetPr(prhs[0]);     // get values
  mwIndex* irs = mxGetIr(prhs[0]);   // get row
  mwIndex* jcs = mxGetJc(prhs[0]);   // get columns
  
  double* Xobs = mxGetPr(prhs[1]);      // Xobs = lambda * s, s = subgradient
  double* rval = mxGetPr(prhs[2]);     // get values    // rval = alpha 

  int MAXITER = mxGetScalar(prhs[3]); 
  double EPS = mxGetScalar(prhs[4]); 
  double L = mxGetScalar(prhs[5]); 
  //mexPrintf("Elements: %f\n",MaxSumSquaredWeights);

  if(L<=0){
	  mexWarnMsgTxt("Lipschitz constant has to be positive");
    return;
  }
  
  double* r2 = mxGetPr(prhs[6]);
//   double* urval = mxGetPr(prhs[7]);
  double* vrval = mxGetPr(prhs[7]);
  double C = mxGetScalar(prhs[8]);
  double normA = mxGetScalar(prhs[9]);
  double *degree = mxGetPr(prhs[10]);
  
//  L = L/(C*C);
  
  //printf("The new variables are set\n");
  double *X;      /* output matrix */
  plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL); /* create the output vector */
  plhs[1] = mxCreateDoubleMatrix(lenrval,1,mxREAL); /* create the output dual variable */
  plhs[2] = mxCreateDoubleMatrix(len,1,mxREAL); /* create the v vector */ 
  plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the output objective value */
  plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); /* create the final iteration value */
  
  //plhs[1]= (mxArray *)prhs[2];
  //plhs[2]= (mxArray *)prhs[3];

  X = mxGetPr(plhs[0]);
  double* Z = mxGetPr(plhs[1]);
  double* V = mxGetPr(plhs[2]);
  double* OutputObj = mxGetPr(plhs[3]);
  double* FinalIter = mxGetPr(plhs[4]);

  int counter=0,i,j,start,mid,end,iter=0;
  double tnew=1; double told=1,alpha,beta,factor,ubeta,vbeta;
  double dummy,normD,normDiff,relativeChange,Fval;
  double* dummyPointer;

  double* D =new double[len];
  double* Dold =new double[len];
  double* Dproj =new double[len];
  for(i=0; i<len; i++) { D[i]=0; Dold[i]=0; Dproj[i]= 0;}

  double* pval    = new double[lenrval];
  double* pvalold = new double[lenrval];
  double* upval = new double[len];
  double* upvalold = new double[len];
  double* vpval = new double[len];
  double* vpvalold = new double[len];
  for(i=0; i<lenrval; i++) { pval[i]=0; }
  for(i=0; i<lenrval; i++) { pvalold[i]=0; }
  for(i=0; i<len; i++) {upval[i] = 0; upvalold[i] = 0; }
  for(i=0; i<len; i++) {vpval[i] = 0; vpvalold[i] = 0; }
  
  //printf("The new variables are initialized\n");
  
  //MaxSumSquaredWeights=max(sum(W.^2,2));
  /*double MaxSumSquaredWeights=0;  
  for(j=0; j<cols; j++) 
  {   
	dummy=0;
	for(i=0; i<jcs[j+1]-jcs[j]; i++) {  dummy+=sr[j]*sr[j]; }
	if(dummy>MaxSumSquaredWeights) { MaxSumSquaredWeights=dummy; }
  }
  MaxSumSquaredWeights=4*MaxSumSquaredWeights;*/
  

  Fval=EPS+1;
//   printf("Fval=%f\t Eps=%f\n", Fval, EPS);
  double primal_obj, max_f, min_f;
  double ratioNormAC = normA/C;
  double ratioWvalC = 0;
    double ratioNormWvalC = 0;
  double scale = 1;
  bool primal_neg = false;
  bool term_cond = false;
//    L = 1;
  double old_Obj = -1000000;  // this assignment not required!
  while(iter<=MAXITER && Fval*C > EPS && !term_cond )// !primal_neg)
//  while(iter<=MAXITER && !term_cond )// !primal_neg)
  {
    // exchange D and Dold
	dummyPointer=D; D=Dold; Dold=dummyPointer;

	//mexPrintf("Exchanged D \n");

	// exchange pval and pvalold
	dummyPointer=pval; pval=pvalold; pvalold=dummyPointer;
    
    // exchange upval and upvalold
	dummyPointer=upval; upval=upvalold; upvalold=dummyPointer;
    
    // exchange vpval and vpvalold
	dummyPointer=vpval; vpval=vpvalold; vpvalold=dummyPointer;

	//mexPrintf("Exchanged Pointer \n");

	// exchange tnew and told
	told=tnew;
    
	// initialize X with zeros 
    for(i=0; i<len; i++) { X[i]=0; }

	//mexPrintf("Initialized X, %i \n",X[0]);
  
    //sval = lambda*wval.*rval;
    //X = -sparse(jx,1,sval);
    dummy=0; counter=0;
    for(j=0; j<cols; j++) 
    {   
	   for(i=0; i<jcs[j+1]-jcs[j]; i++)
	   {  
         //dummy = sr[counter]*rval[counter];
         //dummy = rval[counter]/degree[i];
           dummy = rval[counter];
		 //mexPrintf("Computed dummy, %f\n",dummy);
         X[j] -= dummy;
	     X[irs[counter]] += dummy;
         counter++;
	   }
	} 
    
//     mexPrintf("Aalpha: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", X[i]);
//     mexPrintf("\n");
    
    //mexPrintf("Computed X, %f %f %f %f %f \n",X[0],X[1],X[2],X[3],X[4]);
        
	//D = Xobs+r2-C*u+C*v  -X; 
    //D = Xobs/C +r2/C -u+v  -X/C; 
    //printf("C=%f\tL=%f", C, L);
    
    //D = -Xobs/C -r2/C+u  +X/C; 
    for(i=0; i<len; i++) 
    { 
                //D[i]=Xobs[i]/C+r2[i]/C+vrval[i]-X[i]/C; 
        D[i]=Xobs[i]/C+r2[i]/C+vrval[i]-X[i]/normA; 
                //D[i]=Xobs[i]/C+r2[i]/C+vrval[i]-X[i]*scale; 
//         printf("%f\t", D[i]);
        //D[i] = Xobs[i]-X[i];
//         normD+=D[i]*D[i]; normDiff+=(D[i]-Dold[i])*(D[i]-Dold[i]);
         Dproj[i]=D[i];
    }

//     mexPrintf("b: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", (Xobs[i]+r2[i])/C+vrval[i]);
// //         mexPrintf("%1.15f \t", vrval[i]);
//     printf("\n");   

//         mexPrintf("D: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", D[i]);
// //         mexPrintf("%1.15f \t", vrval[i]);
//     printf("\n");   

    Dproj = ProjKSimplex(Dproj,len,1);
//     Dproj = proj_smplx(Dproj,len);
    
    
//     printf("DProj = ");
    normD = 0;  normDiff=0;
    for(i=0; i<len; i++) 
    { 
//         printf("%f\t", Dproj[i]);
        D[i]=D[i]-Dproj[i];
        normD+=D[i]*D[i]; normDiff+=(D[i]-Dold[i])*(D[i]-Dold[i]);
    }
    
//     mexPrintf("D - Dproj: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", D[i]);
//     printf("\n");   

//     printf("\n");
    
    //printf("normD=%f\n", sqrt(normD));
    //pval = rval + uval.*(D(ix)-D(jx));
	 // pval=pval./max(abs(pval),1);
    counter=0;
	tnew = (1 + sqrt(1+4*told*told))/2;
    factor = (told-1)/tnew;
//     factor = iter/(iter+3);//(told-1)/tnew;
    primal_obj = 0;
    
//     L = 1;
//      mexPrintf("L=%f\n", L);
//         L = BacktrackLineSearch( sr, rval, vrval, D, Xobs, r2, irs, jcs, L, len, lenrval, cols, C, normD, normA, scale, degree);
//         mexPrintf("L=%f\n", L);
//        if( L > 100 )
//          L = L / 10;
       

     //L = sqrt(36);
//      printf("L = %1.15f\n", L);
    
//     mexPrintf("%1.15f\n", normA);
    for(j=0; j<cols; j++) 
    {   
       alpha=D[j];
	   for(i=0; i<jcs[j+1]-jcs[j]; i++)
	   {  
          // beta is x(k+1)    t = 1/MaxSumSquaredWeights
          // alpha = sum_j w_rj
          // D[irs[counter]] = b_r - X_r
	      // update of pval
        		  //beta=rval[counter] + 2*(sr[counter]/normA)*( D[irs[counter]] - alpha)/L;
           beta=rval[counter] + (2/normA)*( D[irs[counter]] - alpha)/L;
           //beta=rval[counter] + 2*(scale/degree[i])*( D[irs[counter]] - alpha)/L;
//           printf("beta=%f\trval[counter]=%f\t", beta, rval[counter]);
          // f = D/norm(D)!
          primal_obj += sr[counter]* fabs( D[irs[counter]] - alpha);
          //beta=rval[counter] + sr[counter]*( D[irs[counter]] )/MaxSumSquaredWeights;
		  // projection onto l_inf-cube 
          ratioNormWvalC = normA*sr[counter]/C;
		  //if(beta>ratioWvalC*scale*degree[i]) beta=ratioWvalC*scale*degree[i];
	      //else if(beta<-ratioWvalC*scale*degree[i]) beta=-ratioWvalC*scale*degree[i];
          if(beta>ratioNormWvalC) beta=ratioNormWvalC;
	      else if(beta<-ratioNormWvalC) beta=-ratioNormWvalC;

		  pval[counter]=beta;
		  // update of rval
           rval[counter] = beta + factor*(beta-pvalold[counter]);
		  counter++;
	   }	  
    }
    
//     mexPrintf("alpha_new: ");
//     for(i=counter-1;i>=counter-5;i--)
//         mexPrintf("%1.15f \t", pval[i]);
//     printf("\n");   
    
    //printf("The alpha part of the variable updated\n"); // Note that we still use old alpha for computing the u, v parts of the next iterate

    //We need urval, upval, upvalold 
    for(i=0;i<len;i++)
    { 
        //upval[i] = urval[i] - 2*C*(C*urval[i] + Dtemp[i]-C*vrval[i])/MaxSumSquaredWeights; 
//        upval[i] = urval[i] - 2*(D[i])/L; 
         vpval[i] = vrval[i] - 2*(D[i])/L; 
    }
//     mexPrintf("v: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", vpval[i]);
//     printf("\n");   
    
//     double *samp;
//     double samp1[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
//     
//     samp = ProjKSimplex(samp1, 5, 1);
//     for(i=0;i<5;i++){
//        printf("%f\t ", samp[i]);
//     }
//     printf("%\n");
    
	//projection onto the simplex
    //upval = ProjKSimplex(upval, len, 1);
    vpval = ProjKSimplex(vpval, len, 1);
//     vpval = proj_smplx(vpval, len);
//     mexPrintf("vp: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", vpval[i]);
//     printf("\n");   
    

	// update y(k+1)
    for(i=0;i<len;i++)
    {
   //     urval[i] = upval[i] + factor*(upval[i]-upvalold[i]);
        vrval[i] = vpval[i] + factor*(vpval[i]-vpvalold[i]);
    }
    
//     mexPrintf("vr: ");
//     for(i=0;i<5;i++)
//         mexPrintf("%1.15f \t", vrval[i]);
//     printf("\n");   

    //tkp1=(1+sqrt(1+4*tk^2))/2;
    //rval = pval + (tk-1)/(tkp1)*(pval-pvalold);
	/*tnew = (1 + sqrt(1+4*told*told))/2;
    for(j=0; j<jcs[len]; j++) 
    {
	  alpha=pval[j];
	  if(alpha>1) { pval[j]=1; alpha=1;  }
	  else if(alpha<-1) { pval[j]=-1; alpha=-1;}
      rval[j] = alpha + (told-1)/tnew*(alpha-pvalold[j]);
    }*/
	//mexPrintf("Comp pval: %f %f %f\n",pval[0],pval[1],pval[len-1]);
	//mexPrintf("Comp rval: %f %f %f\n",rval[0],rval[1],rval[len-1]);
  
    relativeChange = sqrt(normDiff/normD);

    Fval = sqrt(normD);
    
        
    for(i=0;i<len;i++)
        primal_obj += -D[i]*Xobs[i];
    
    for(i=0;i<len;i++)
        primal_obj += -D[i]*r2[i];
    
    // Find max and min of f
//     printf("f = ");
    //min_f = DBL_MAX; max_f = -min_f;//min(double);
    min_f = D[0]; 
    max_f = D[0];
    for(i=1; i<len; i++)
    {
        
//         printf("%f\t", D[i]);
        if(D[i] < min_f)
            min_f = D[i];
        if(D[i] > max_f)
            max_f = D[i];
    }
//     printf("\nmax_f = %f\t min_f=%f\n", max_f, min_f);
    
    primal_obj += C*max_f;
    primal_obj += -C*min_f;
    primal_obj/= sqrt(normD);
    
//     if( -Fval*C < old_Obj)
//     {
//         mexPrintf("dual objective decreases... bug in fista code: old_Ojb = %1.16f \t Obj = %1.16f\n", old_Obj, -Fval*C);
//     }
            

//     mexPrintf("normD=%f\n", normD);
    
//     if( (primal_obj+Fval*C) < 0 ) // Check the duality gap
//     {
// //         term_cond = true;
//         mexPrintf("bug in fista code???? is it infeasible dual??\n");
// //         primal_obj = DBL_MAX;
// //         Fval = DBL_MAX;
//     }
    
    double EPS1 = 1e-3;
    if (L >= 1/EPS || fabs(primal_obj+Fval*C)/(Fval*C) < EPS)
    {
        term_cond = true;
//         mexPrintf("stopping the inner iterations since the duality gap = %1.16f is small (%1.15f) or step size L = %1.16f is too large \n", (primal_obj+Fval*C)/(Fval*C), EPS1, L);
    }
    
    
//     if L > 4

//     if(primal_obj < 0)
//         primal_neg = true;


    //mexPrintf("Iteration: %i, Fval: %1.15f, RelativeChange %1.15f\n",iter,Fval,relativeChange);
//  	if((iter>0 &&iter<10) || (iter>10 &&iter % 10==0))
// 	 mexPrintf("Iteration: %i, Fval: %1.15f, RelativeChange %1.15f\n",iter,Fval,relativeChange);
//     if(iter<10 || iter % 100==0) //|| iter == MAXITER-1)
 	 //mexPrintf("Iteration: %i, Fval: %1.15f, RelativeChange %1.15f\n",iter,Fval,relativeChange);
//     mexPrintf("Iteration: %i, Dual: %1.15f, Primal: %1.15f, RelativeGap %1.15f\tL=%1.15f\n",iter,-Fval*C,primal_obj, (primal_obj+Fval*C)/(Fval*C),L);
//      mexPrintf("Iteration: %i, Dual: %1.15f, Primal: %1.15f, RelativeGap %1.15f\tL: %1.15f\n",iter,-Fval*C,primal_obj, (primal_obj+Fval*C)/(Fval*C), L);
    
    
//     double mini = 2*fabs(old_Obj + Fval*C);
// //     mexPrintf("mini=%f\n", mini);
//     if (mini > 1)
//         mini = 1;
//    
//     L = mini/(8*normD);
//     L = 1/L;

    old_Obj = -Fval*C;
	iter++;
 }
// mexPrintf("FINAL: Iterations %i, Fval: %1.15f,x RelativeChange %1.15f\n",iter,Fval,relativeChange);
  
//   mexPrintf("Iteration: %i, Dual: %1.15f, Primal: %1.15f, RelativeGap %1.15f\tL=%1.15f\n",iter,-Fval*C,primal_obj, (primal_obj+Fval*C)/(Fval*C),L);
//   mexPrintf("Iteration: %i, Dual: %1.15f, Primal: %1.15f, RelativeGap %1.15f\tL: %1.15f\n",iter,-Fval*C,primal_obj, (primal_obj+Fval*C)/(Fval*C), L);
 for(i=0; i<len; i++) { X[i]=D[i];}
 for(i=0; i<len; i++) { V[i]=vrval[i];}
 for(i=0; i<lenrval; i++) { Z[i]=rval[i];}
 OutputObj[0]=Fval*C;
 FinalIter[0]=iter;
   
 //delete D; delete Dold; delete pvalold; delete pval;
 // Delete all memory created using new!
 delete D; delete Dold; delete Dproj;
  delete pvalold;  delete pval; 
  delete upval; delete upvalold; 
  delete vpval; delete vpvalold;
}


double BacktrackLineSearch( double *sr, double *alpha, double *v, double *D, double* Xobs, double* r2, mwIndex *irs, mwIndex*jcs, double L, int len, int m, int cols, double C, double old_obj, double normA, double scale, double *degree)
{

    double mu = 2.1;    
//     double L = L_old;
    bool flag = true;
    
    double* p_alpha = new double[m];
    double* pv = new double[len];
    double* Dtemp = new double[len];
    double* Dprojtemp = new double[len];
    double* X = new double[len];
    
    double* diff_alpha = new double[m];
    double* diff_v = new double[len];
    
    int counter = 0;
    int i, j;
    double temp=0, beta, dummy, lhs, rhs=0, qterm;
    double normD = 0;
    double iter = 0;
    double ratioNormAC = normA/C;
    double ratioWvalC;
    double ratioNormWvalC;
//     double scale = 10;
    
//     printf("D = ");
//     for(i=0;i<len;i++){
//         printf("%f\t", D[i]);
//         if(i%8 == 0)
//             printf("\n");
//     }
// 	printf("\nv = ");
//     for(i=0;i<len;i++){
//         printf("%f\t", v[i]);
//         if(i%8 == 0)
//             printf("\n");
//     }
//     printf("\n\n\n");
    
//     counter = 0;
//     for(j=0;j<cols;j++)
//         for(i=0;i<jcs[j+1]-jcs[j];i++)
//             printf("%f\t", sr[counter++]);
//     printf("\n");
//     printf("old_ojb=%f\n",old_obj);
//     L = 10;
    while (flag)
    {
//         
// %          alpha_new = alpha + (1/L) * descent_direction;
//         %p_alpha = Projection( alpha + (1/L)*descent_direction );
//         
//         p_alpha = alpha - (1/L)*grad_alpha;
//         p_alpha(p_alpha>1) = 1;
//         p_alpha(p_alpha<-1) = -1;
//         %p_u = u - (1/L)*grad_u;
//         p_v = v - (1/L)*grad_v;
//         p_v = projectSimplex(p_v);
        
        for(i=0; i<len; i++) { X[i]=0; }
        counter = 0;
        double beta_temp;
        for(j=0; j<cols; j++) 
        {   
           temp=D[j];
           for(i=0; i<jcs[j+1]-jcs[j]; i++)
    	   {  
              beta=alpha[counter] + 2*(sr[counter]/normA)*( D[irs[counter]] - temp)/L; 
               //beta=alpha[counter] + 2*(scale/degree[i])*( D[irs[counter]] - temp)/L; 
//               beta_temp = beta;
//               printf("%f\t", beta);
              // projection onto l_inf-cube 
               ratioNormWvalC = normA*sr[counter] / C;
              if(beta>ratioNormWvalC) beta=ratioNormWvalC;
              else if(beta<-ratioNormWvalC) beta=-ratioNormWvalC;
               //if(beta>ratioWvalC*scale*degree[i]) beta=ratioWvalC*scale*degree[i];
              //else if(beta<-ratioWvalC*scale*degree[i]) beta=-ratioWvalC*scale*degree[i];
              p_alpha[counter]=beta;
              diff_alpha[counter] = p_alpha[counter] - alpha[counter];
//               printf("diff = %f\t", diff_alpha[counter]);
              counter++;
           }	              
        }
//         printf("rationwvalC = %1.16f\n", ratioWvalC);
        
        for(i=0;i<len;i++)
        { 
             pv[i] = v[i] - 2*(D[i])/L; 
        }
        pv = ProjKSimplex(pv, len, 1);
//         pv = proj_smplx(pv, len);
        for(i=0;i<len;i++){
            diff_v[i] = pv[i] - v[i];
//             printf("%f\t", diff_v[i]);
            
        }

        //         
//         Apalpha = sparse(ix,1,wval.*p_alpha, n, 1) - sparse(jx,1,wval.*p_alpha, n,1);    
//         pD = s/C+r2/C+p_v - Apalpha/C;
//         pD = pD - projectSimplex(pD);

        dummy=0; counter=0;
        for(j=0; j<cols; j++) 
        {   
           for(i=0; i<jcs[j+1]-jcs[j]; i++)
           {  
             //dummy = sr[counter]*p_alpha[counter];
               dummy = p_alpha[counter];
              // dummy = p_alpha[counter]/degree[i];
             //mexPrintf("Computed dummy, %f\n",dummy);
             X[j] -= dummy;
             X[irs[counter]] += dummy;
//              printf("%f\t", dummy);
             counter++;
           }
        } 
        
        for(i=0; i<len; i++) 
        { 
            Dtemp[i]=Xobs[i]/C+r2[i]/C+pv[i]-X[i]/normA; 
            //Dtemp[i]=Xobs[i]/C+r2[i]/C+pv[i]-X[i]*scale; 
//             sprintf("%f\t", X[i]);
            Dprojtemp[i]=Dtemp[i];
        }

        Dprojtemp = ProjKSimplex(Dprojtemp,len,1);
//         Dprojtemp = proj_smplx(Dprojtemp,len);

        normD = 0;
        for(i=0; i<len; i++) 
        { 
            Dtemp[i]=Dtemp[i]-Dprojtemp[i];
            //printf("%f\t", Dtemp[i]);
              normD+=Dtemp[i]*Dtemp[i]; 
//             normDiff+=(D[i]-Dold[i])*(D[i]-Dold[i]);
        }
        
        //         
//         lhs = norm(pD)^2;
//         %pvec = [p_alpha-alpha; p_u- u; p_v-v];
//         pvec = [p_alpha-alpha; p_v-v];
//         
//         %rhs = old_obj + (pvec)'* [grad_alpha; grad_u; grad_v] + (L/2)* norm(pvec,2)^2;
//         rhs = old_obj + (pvec)'* [grad_alpha; grad_v] + (L/2)* norm(pvec,2)^2;
// end
        
        lhs = normD;
        rhs = old_obj;
//         printf("old_obj = %1.16f\n", old_obj);
        // inner product: linear term
        //for(i=0;i<counter-1;i++)
        counter = 0;
        qterm = 0;
        double ft1 = 0, ft2 = 0;
        for(j=0; j<cols; j++) 
        {   
           temp=D[j];
           for(i=0; i<jcs[j+1]-jcs[j]; i++)
    	   {  
               rhs += diff_alpha[counter] * 2*(-sr[counter]/normA)*( D[irs[counter]] - temp);
               //rhs += diff_alpha[counter] * -2*(scale/degree[i])*( D[irs[counter]] - temp);
               ft1 += diff_alpha[counter] * 2*(-sr[counter]/normA)*( D[irs[counter]] - temp);
               //ft1 += diff_alpha[counter] * -2*(scale/degree[i])*( D[irs[counter]] - temp);
               qterm += diff_alpha[counter]*diff_alpha[counter];
               counter++;
           }
        }
        
        for(i=0;i<len;i++)
        {
             rhs += diff_v[i] * 2*(D[i]); 
             ft1 += diff_v[i] * 2*(D[i]); 
             ft2 += diff_v[i]*diff_v[i]; 
            qterm += diff_v[i]*diff_v[i]; 
        }
        qterm *= L/2;
         rhs += qterm;
//           printf("%1.16f %1.16f %1.16f\n", ft1, ft2, L);
        
        //         if lhs > rhs
//             L = beta * L;
// %            p_alpha = Projection( alpha + (1/L)*descent_direction );
//         else
//             flag = false;
//         end
//     end
// 

//             printf( "L = %1.16f \t lhs = %1.16f \t rhs = %1.16f\n", L, lhs, rhs);
        if (lhs > rhs)
            L = mu*L;
        else
            flag = false;    
        
        if (iter > 200000){
            mexPrintf("\n\n\nBug in the line search\n\n\n\n");
            flag = false;
        }
        
        iter ++; 

    }
   
    delete p_alpha; delete pv; delete Dtemp; delete Dprojtemp; 
    delete X; delete diff_alpha; delete diff_v;

    //printf("L = %f\n", L);
//     printf("L=%1.15f\n",L);
    return L;
}


double* ProjKSimplex(double* f,int len, int K)
{
  int i;
  //double* xnew = new double[len];
  int* IX = new int[len];
  double sumxnew = 0; double sumxnew2;
  int counter=len; int counter2=0;
  for(i=0; i<len; i++)
   { IX[i]=i; sumxnew+=f[i];}
  int iteration=0;
  while(iteration<10000)
  {
    counter2=0; sumxnew2=0;
    for(i=0;i<counter; i++)
	{
	  f[IX[i]] = f[IX[i]] - (sumxnew-K)/counter;
	  if(f[IX[i]]>0)
	   { IX[counter2]=IX[i]; counter2++; sumxnew2+=f[IX[i]]; }
	  else
	   f[IX[i]]=0;
	}
	if(counter2==counter) break;
	counter=counter2; sumxnew=sumxnew2;
	iteration++;
  }
  if(iteration==10000)
	  mexPrintf("Bug in Projection: %i, Fval: %i\n",counter,counter2);
  delete IX;
  return f;
  /*x = f;  IX=1:num;  
  while(cond)
   xnew = x;
   xnew(IX) = x(IX) - (sum(x)-k)/length(IX);
   if(sum(xnew>=0)==num)
    Projf=xnew; break;
   else
    %ixx = find(xnew<0);
    %IX = setdiff(IX,ixx);
    IX = find(xnew>0);
    x = max(xnew,0);
   end
  end*/
}


// typedef double ElementType;
//      
//      
// void Swap( ElementType *Lhs, ElementType *Rhs )
// {
//             ElementType Tmp = *Lhs;
//             *Lhs = *Rhs;
//             *Rhs = Tmp;
// }
// ElementType Median3( ElementType A[ ], int Left, int Right )
// {
//             int Center = ( Left + Right ) / 2;
// 
//             if( A[ Left ] > A[ Center ] )
//                 Swap( &A[ Left ], &A[ Center ] );
//             if( A[ Left ] > A[ Right ] )
//                 Swap( &A[ Left ], &A[ Right ] );
//             if( A[ Center ] > A[ Right ] )
//                 Swap( &A[ Center ], &A[ Right ] );
// 
//             /* Invariant: A[ Left ] <= A[ Center ] <= A[ Right ] */
// 
//             Swap( &A[ Center ], &A[ Right - 1 ] );  /* Hide pivot */
//             return A[ Right - 1 ];                /* Return pivot */
// }
// void InsertionSort( ElementType A[ ], int N )
// {
//             int j, P;
//             ElementType Tmp;
// 
// /* 1*/      for( P = 1; P < N; P++ )
//             {
// /* 2*/          Tmp = A[ P ];
// /* 3*/          for( j = P; j > 0 && A[ j - 1 ] > Tmp; j-- )
// /* 4*/              A[ j ] = A[ j - 1 ];
// /* 5*/          A[ j ] = Tmp;
//             }
// }
//         
// #define Cutoff ( 3 )
// 
// void Qsort( ElementType A[ ], int Left, int Right )
// {
//             int i, j;
//             ElementType Pivot;
//             if( Left + Cutoff <= Right )
//             {
// 	        Pivot = Median3( A, Left, Right );
//                 i = Left; j = Right - 1;
// 	        for( ; ; )
//                 {
// 	            while( A[ ++i ] < Pivot ){ }
// 	            while( A[ --j ] > Pivot ){ }
//                     if( i < j )
//                        Swap( &A[ i ], &A[ j ] );
//                     else
//                        break;
//                 }
//                 Swap( &A[ i ], &A[ Right - 1 ] );  /* Restore pivot */
// 
//                 Qsort( A, Left, i - 1 );
//                 Qsort( A, i + 1, Right );
//             }
//             else  /* Do an insertion sort on the subarray */
// 	        InsertionSort( A + Left, Right - Left + 1 );
// }
//       
//           
// void Quicksort( ElementType A[ ], int N )
// {
//             Qsort( A, 0, N - 1 );
// }
// 
// double* proj_smplx( double* y, double m )
// {
//     int numDims, d, j;
//     //int *dims;
//     // double *y, *s, *x;
//     double *x, *s;
//     double sumResult = -1, tmpValue, tmax; 
//     bool bget = false;
//  
//     //printf("inside the projection function\n");
// //     dims = mxGetDimensions(prhs[0]);
// //     numDims = mxGetNumberOfDimensions(prhs[0]);  
// 
//     /*m = dims[0]; n=dims[1];*/
// //     m = 0;
// //     for (d=0; d<numDims; d++) {
// //        m= (m > dims[d])? m:dims[d];
// //     }
// 
//     //y  = mxGetPr(prhs[0]); 
//      
//      
//     /*  set the output pointer to the output matrix */
//     //plhs[0] = mxCreateDoubleMatrix(m,1,mxREAL);
//     x = (double*) calloc (m,sizeof(double));
//     
//     /* s = sort(y,'ascend'); */
//     s = (double*) calloc (m,sizeof(double));
//     for(j = 0; j < m; j++ ){
//     	s[j] = y[j]; 
//     }
//     Quicksort(s,m);
//     
//     //x = mxGetPr(plhs[0]);
// 
//     /* if t is not less than s[0] */
//     for(j = m-1; j >= 1; j--){    	
//     	sumResult = sumResult + s[j];
//     	tmax = sumResult/(m-j);
// 	if(tmax >= s[j-1]){
// 		bget = true;
// 		break;
// 	}
//     }   
// 
//     /* if t is less than s[0] */
//     if(!bget){
// 	sumResult = sumResult + s[0];
// 	tmax = sumResult/m;
//     }
// 
//     /* x = max(y-tmax, 0); */
//     for(j = 0; j <= m-1; j++){
// 	tmpValue = y[j] - tmax;
// 	x[j] = (tmpValue > 0)? tmpValue:0;
//     }
//     
//     //printf("projection onto simplex completed\n");
//     return x;
// }
// 
