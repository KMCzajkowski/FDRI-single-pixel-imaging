    #include "mex.h"

void
mexFunction (int nlhs, mxArray* plhs[],
             int nrhs, const mxArray* prhs[])
{
  mwSize n;
  mwIndex i;
  double *vri, *vro, *vro2;
    
  double binval,oldval,newval,var,del;

  if (nlhs != 2 || ! mxIsDouble (prhs[0]))
    mexErrMsgTxt ("BBB");

  n = mxGetNumberOfElements (prhs[0]);
  plhs[0] = mxCreateNumericArray (mxGetNumberOfDimensions (prhs[0]),
                                  mxGetDimensions (prhs[0]),
                                  mxGetClassID (prhs[0]),
                                  mxIsComplex (prhs[0]));
  plhs[1] = mxCreateNumericArray (mxGetNumberOfDimensions (prhs[0]),
                                  mxGetDimensions (prhs[0]),
                                  mxGetClassID (prhs[0]),
                                  mxIsComplex (prhs[0]));
  vri = mxGetPr (prhs[0]);
  vro = mxGetPr (plhs[0]);
  vro2 = mxGetPr (plhs[1]);
  
  vro[0]=0;
  vro[n]=0;
  vro2[0]=0;
  vro2[n]=0;
  newval=vri[0];
  var=7.0/16.0;
  for(i=1;i<(n-1);i++){
        oldval=newval;
          binval=0;
          if(oldval>0.5){
            binval=1;
          }
        vro[i]=binval;
        del=oldval-binval;

        newval= vri[i+1]+del*var;
        vro2[i]=del;
        }
}
