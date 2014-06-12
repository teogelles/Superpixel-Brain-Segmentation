#include <math.h>
#include "mex.h"
#include "UGM_common.h"
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Variables */
  int n, s,e,e2,n1,n2,neigh,Vind,Vind2,s1,s2,
    nNodes, nEdges, maxState, dims[3],
    iter,maxIter,nNbrs,
    *edgeEnds, *nStates, *V, *E,*y;
    
  double *nodePot, *edgePot, *nodeBel, *edgeBel, *logZ,
    z,energy1,energy2,entropy1,entropy2,*prodMsgs,*oldMsgs,*newMsgs;
    
  /* Input */
    
  nodePot = mxGetPr(prhs[0]);
  edgePot = mxGetPr(prhs[1]);
  edgeEnds = (int*)mxGetPr(prhs[2]);
  nStates = (int*)mxGetPr(prhs[3]);
  V = (int*)mxGetPr(prhs[4]);
  E = (int*)mxGetPr(prhs[5]);
  maxIter = ((int*)mxGetPr(prhs[6]))[0];
    
  if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32"))
    mexErrMsgTxt("edgeEnds, nStates, V, E, maxIter must be int32");
	
  /* Compute Sizes */
    
  nNodes = mxGetDimensions(prhs[0])[0];
  maxState = mxGetDimensions(prhs[0])[1];
  nEdges = mxGetDimensions(prhs[2])[0];
    
  /* Output */
  plhs[0] = mxCreateDoubleMatrix(nNodes,maxState,mxREAL);
  dims[0] = maxState;
  dims[1] = maxState;
  dims[2] = nEdges;
  plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
  nodeBel = mxGetPr(plhs[0]);
  edgeBel = mxGetPr(plhs[1]);
  logZ = mxGetPr(plhs[2]);
    
  prodMsgs = mxCalloc(maxState*nNodes,sizeof(double));
  oldMsgs = mxCalloc(maxState*nEdges*2,sizeof(double));
  newMsgs = mxCalloc(maxState*nEdges*2,sizeof(double));

  /* Initialize */
  for(e = 0; e < nEdges; e++)
    {
      n1 = edgeEnds[e]-1;
      n2 = edgeEnds[e+nEdges]-1;
      for(s = 0; s < nStates[n2]; s++)
	newMsgs[s+maxState*e] = 1./nStates[n2];
      for(s = 0; s < nStates[n1]; s++)
	newMsgs[s+maxState*(e+nEdges)] = 1./nStates[n1];
    }


  (void) omp_set_dynamic(0);
  (void) omp_set_num_threads((omp_get_num_procs()));
    

  for(iter = 0; iter < maxIter; iter++)
    {
#pragma omp parallel for shared(nodePot,edgePot,V,E,edgeEnds,nStates,newMsgs,nEdges,maxState,nNodes) default(none) private(Vind,Vind2,s,s1,s2,n,n1,n2,e,e2,z)
      for(n=0;n<nNodes;n++)
        { 
	  double *tmp = calloc(maxState,sizeof(double));
	  /* Update Messages */
	  for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
            {                
	      e = E[Vind]-1;
	      n1 = edgeEnds[e]-1;
	      n2 = edgeEnds[e+nEdges]-1;
                
	      /* First part of message is nodePot*/
	      for(s = 0; s < nStates[n]; s++)
                {
		  tmp[s] = nodePot[n + nNodes*s];
                }
                
	      /* Multiply by messages from neighbors except j */
	      for(Vind2 = V[n]-1; Vind2 < V[n+1]-1; Vind2++)
                {
		  e2 = E[Vind2]-1;
		  if (e != e2)
                    {
		      if (n == edgeEnds[e2+nEdges]-1)
                        {
			  for(s = 0; s < nStates[n]; s++)
                            {
			      tmp[s] *= newMsgs[s+maxState*e2];
                            }
                        }
		      else
                        {
			  for(s = 0; s < nStates[n]; s++)
                            {
			      tmp[s] *= newMsgs[s+maxState*(e2+nEdges)];
                            }
                        }
                    }
                }
                
                
	      /* Now multiply by edge potential to get new message */
                
	      if (n == n2)
                {
                    
		  for(s1 = 0; s1 < nStates[n1]; s1++)
                    {
		      newMsgs[s1+maxState*(e+nEdges)] = 0.0;
		      for(s2 = 0; s2 < nStates[n2]; s2++)
                        {
			  newMsgs[s1+maxState*(e+nEdges)] += tmp[s2]*edgePot[s1+maxState*(s2+maxState*e)];
                        }
                        
                    }
                    
		  /* Normalize */
		  z = 0.0;
		  for(s = 0; s < nStates[n1]; s++)
		    z += newMsgs[s+maxState*(e+nEdges)];
		  for(s = 0; s < nStates[n1]; s++)
		    newMsgs[s+maxState*(e+nEdges)] /= z;
                }
	      else
                {
		  for(s2 = 0; s2 < nStates[n2]; s2++)
                    {
		      newMsgs[s2+maxState*e] = 0.0;
		      for(s1 = 0; s1 < nStates[n1]; s1++)
                        {
			  newMsgs[s2+maxState*e] += tmp[s1]*edgePot[s1+maxState*(s2+maxState*e)];
                        }
                        
                    }
                    
		  /* Normalize */
		  z = 0.0;
		  for(s = 0; s < nStates[n2]; s++)
		    z += newMsgs[s+maxState*e];
		  for(s = 0; s < nStates[n2]; s++)
		    newMsgs[s+maxState*e] /= z;
                }
                
            }            
	  free(tmp);               
            
            
        }
        
      /* oldMsgs = newMsgs */
      z = 0;
      for(s=0;s<maxState;s++)
        {
#pragma omp parallel for default(shared) private(e) reduction(+:z)
	  for(e=0;e<nEdges*2;e++)
            {
	      z += absDif(newMsgs[s+maxState*e],oldMsgs[s+maxState*e]);
	      oldMsgs[s+maxState*e] = newMsgs[s+maxState*e];
            }
        }
        
      /* if sum(abs(newMsgs(:)-oldMsgs(:))) < 1e-4; break; */
      if(z < 1e-4)
        {
	  break;
        }
    }
    
  /*if(iter == maxIter)
    {
    printf("LBP reached maxIter of %d iterations\n",maxIter);
    }
    printf("Stopped after %d iterations\n",iter); */
    
  /* compute nodeBel */
#pragma omp parallel for default(none) shared(edgeEnds,nEdges,prodMsgs,maxState,nodePot,V,E,newMsgs,nodeBel,nStates,nNodes) private(s,n,Vind,e,n1,n2,z)
  for(n = 0; n < nNodes; n++)
    {
      for(s = 0; s < nStates[n]; s++)
	prodMsgs[s+maxState*n] = nodePot[n+nNodes*s];
        
      for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
        {
	  e = E[Vind]-1;
	  n1 = edgeEnds[e]-1;
	  n2 = edgeEnds[e+nEdges]-1;
            
	  if (n == n2)
            {
	      for(s = 0; s < nStates[n]; s++)
                {
		  prodMsgs[s+maxState*n] *= newMsgs[s+maxState*e];
                }
            }
	  else
            {
	      for(s = 0; s < nStates[n]; s++)
                {
		  prodMsgs[s+maxState*n] *= newMsgs[s+maxState*(e+nEdges)];
                }
            }
        }
        
      z = 0;
      for(s = 0; s < nStates[n]; s++)
        {
	  nodeBel[n + nNodes*s] = prodMsgs[s+maxState*n];
	  z = z + nodeBel[n+nNodes*s];
        }
      for(s = 0; s < nStates[n]; s++)
	nodeBel[n + nNodes*s] /= z;
    }
    
#pragma omp parallel for default(none) shared(maxState,edgePot,V,E,newMsgs,edgeEnds,nodeBel,nNodes,edgeBel,nStates,nEdges) private(s,s1,s2,n,Vind,e,n1,n2,z)
    
  /* Compute edgeBel */
  for(e = 0; e < nEdges; e++)
    {
      n1 = edgeEnds[e]-1;
      n2 = edgeEnds[e+nEdges]-1;
      z = 0;
      for(s1 = 0; s1 < nStates[n1]; s1++)
        {
	  for(s2 = 0; s2 < nStates[n2]; s2++)
            {
	      edgeBel[s1+maxState*(s2+maxState*e)] = nodeBel[n1+nNodes*s1]/newMsgs[s1+maxState*(e+nEdges)];
	      edgeBel[s1+maxState*(s2+maxState*e)] *= nodeBel[n2+nNodes*s2]/newMsgs[s2+maxState*e];
	      edgeBel[s1+maxState*(s2+maxState*e)] *= edgePot[s1+maxState*(s2+maxState*e)];
	      z += edgeBel[s1+maxState*(s2+maxState*e)];
            }
        }
      for(s1 = 0; s1 < nStates[n1]; s1++)
        {
	  for(s2 = 0; s2 < nStates[n2]; s2++)
	    edgeBel[s1+maxState*(s2+maxState*e)] /= z;
        }
    }

  /* Compute Bethe Free Energy */
  energy1 = 0;
  energy2 = 0;
  entropy1 = 0;
  entropy2 = 0;
#pragma omp parallel for shared(nNodes,nStates,nodeBel,nodePot,V) default(none) private(nNbrs,s,n) reduction(+:entropy1) reduction(-:energy1)
  for(n = 0; n < nNodes; n++)
    {
      nNbrs = V[n+1]-V[n];
      for(s = 0; s < nStates[n]; s++)
        {
	  if(nodeBel[n+nNodes*s] > 1e-10)
	    entropy1 += (nNbrs-1)*nodeBel[n+nNodes*s]*log(nodeBel[n+nNodes*s]);
            
	  energy1 -= nodeBel[n+nNodes*s]*log(nodePot[n+nNodes*s]);
        }
    }

#pragma omp parallel for shared(nEdges,nStates,edgeEnds,edgeBel,edgePot,maxState) default(none) private(e,n1,n2,s1,s2) reduction(-:energy2,entropy2)
  for(e = 0; e < nEdges; e++)
    {
      n1 = edgeEnds[e]-1;
      n2 = edgeEnds[e+nEdges]-1;
        
      for(s1 = 0; s1 < nStates[n1];s1++)
        {
	  for(s2 = 0; s2 < nStates[n2]; s2++)
            {
	      if(edgeBel[s1+maxState*(s2+maxState*e)] > 1e-10)
                {
		  entropy2 -= edgeBel[s1+maxState*(s2+maxState*e)]*log(edgeBel[s1+maxState*(s2+maxState*e)]);
                }
	      energy2 -= edgeBel[s1+maxState*(s2+maxState*e)]*log(edgePot[s1+maxState*(s2+maxState*e)]);
            }
        }
    }
  logZ[0] = -energy1-energy2+entropy1+entropy2;
    
  /* Free memory */
  mxFree(prodMsgs);
  mxFree(oldMsgs);
  mxFree(newMsgs);
}
