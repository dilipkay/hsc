#include "mex.h"
#include<iostream>
#include <math.h>
#include <time.h>
#include<map>
#include "mxGetPropertyPtr.h"
#include "multilevel.h"
using namespace std;

/*
  Given a matrix A and division of the variables into
  coarse and fine, determine a prolongation matrix based
  on Gaussian elimination of the fine variables.
 */

// Assumed that A is sparse and C is a vector
void hsc_calc_prolong(mxArray *A, mxArray *C, mxArray *P, mxArray *invD) 
{
	mwIndex *A_irow = mxGetIr(A);
	mwIndex *A_jcol = mxGetJc(A);
	double *A_pr = (double *) mxGetPr(A);
	const mwIndex *dim = mxGetDimensions(A);
	int n = dim[0];
	const mwIndex *dim1 = mxGetDimensions(P);
	int m = dim1[1];
	unsigned *C_pr = (unsigned *) mxGetPr(C);
	mwIndex *P_irow = mxGetIr(P);
	mwIndex *P_jcol = mxGetJc(P);
	double *P_pr = (double *) mxGetPr(P);
	double *invD_pr = (double *) mxGetPr(invD);
	clock_t tt1, tt2, st;

	tt1 = clock();
	int nz = 0, cc = 0;
	for (int i=0; i<n; i++) {
		if(C_pr[i] == 1) { // coarse point
			P_jcol[cc] = nz;
      
			for (int j=A_jcol[i]; j<A_jcol[i+1]; j++) {
			// only consider non-zero neighbors
			//if (A_pr[j] >= 0) 
			// continue; 
			int row = A_irow[j];

			P_irow[nz] = row;

			if (C_pr[row] == 0)  
				P_pr[nz++] = -A_pr[j]*invD_pr[row];
			else if (row == i){
				P_pr[nz++] = 1;
			}
		}
		cc++;
		}
	}
	tt2 = clock();
  
	if (cc != m)
		fprintf(stderr, "This is bad news in hsc_calc_prolong that cc != m cc=%d m=%d\n", cc, m);

	P_jcol[m] = nz;
	fprintf(stderr, "Computing prolongation matrix %f s\n", double(tt2-tt1)/CLOCKS_PER_SEC);
}
