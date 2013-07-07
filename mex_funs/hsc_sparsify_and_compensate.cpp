#include "mex.h"
#include<iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mxGetPropertyPtr.h"
#include "multilevel.h"
using namespace std;

/*
  Given an M-matrix A, use energy-based sparsification techniques
  and greedy front-growing to create select the coarse level variables
  and create a restriction/prolongation matrix.
 */

#define A_const prhs[0]
#define X prhs[1] // X coordinates of variables 
#define Y prhs[2] // Y coordinates of variables 
#define Z prhs[3] // Z coordinates of variables (0 for 3D grid)
#define Level prhs[4] // current Level (starting at 1 for finest)
#define Ge_thresh prhs[5] // threshold that selects geometric edge or weakest edge to cut in each triangle 
#define Geom_info prhs[6] // flag to determine if geometric information can be used
#define Multi_triangle prhs[7] // flag to determine if multi-triangle compensation to be done during sparsification
#define P plhs[1] // output prolongation matrix
#define C plhs[2] // output coarse indices
#define invD plhs[3] // inverse of A's diagonal 
#define tsp plhs[4] // sparsification time

extern "C" bool mxUnshareArray(mxArray *array_ptr, bool noDeepCopy);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	mxArray *A = const_cast<mxArray *>(A_const);
	mxUnshareArray(A, true);
	plhs[0] = A;

	const mwIndex *dim = mxGetDimensions(A);

	int n = dim[0];
	clock_t st, et;
  
	st = clock();
 
	mwIndex *A_irow = mxGetIr(A);
	mwIndex *A_jcol = mxGetJc(A);
	double *A_pr = (double *) mxGetPr(A);
	mwSize dims[1]; dims[0] = n;
  
	int level = (int) (*mxGetPr(Level));
	double ge_thresh = *mxGetPr(Ge_thresh);
	unsigned geom_info = (unsigned) *mxGetPr(Geom_info);
	unsigned multi_triangle = (unsigned) *mxGetPr(Multi_triangle);


	mwSize dims1[1]; 
	dims1[0] = 1;
	tsp = mxCreateNumericArray(1, dims1, mxDOUBLE_CLASS, mxREAL);
	double *tsp_pr = (double *) mxGetPr(tsp);
  
	int rows, cols;

	double *X_pr, *Y_pr, *Z_pr;
	if (geom_info == 0) {
		X_pr = NULL;
		Y_pr = NULL;
		Z_pr = NULL;
	} else if (geom_info == 1) {
    	X_pr = (double *) mxGetPr(X); 
    	Y_pr = (double *) mxGetPr(Y);
    	Z_pr = NULL;
	} else if (geom_info == 2) {
    	X_pr = (double *) mxGetPr(X); 
    	Y_pr = (double *) mxGetPr(Y);
    	Z_pr = (double *) mxGetPr(Z);
	}

	invD = mxCreateNumericArray(1, dims, mxDOUBLE_CLASS, mxREAL);
	double *invD_pr = (double *) mxGetPr(invD);

	int k = 0;

	// select coarse variables
	C = mxCreateNumericArray(1, dims, mxUINT32_CLASS, mxREAL);
	int num_coarse;

  	// sparsify
	num_coarse = hsc_sparsify(A, geom_info, X_pr, Y_pr, Z_pr, ge_thresh, C, level - 1, tsp_pr, invD_pr, multi_triangle);

	fprintf(stderr, "Level %d coarse variables = %f%%\n\n", level, 100.0*float(num_coarse)/float(n));

	// create prolongation matrix
	P = mxCreateSparse(n, num_coarse, mxGetNzmax(A), mxREAL);

	int num_fine = n - num_coarse;
	// Calculate the prolongation matrix
	hsc_calc_prolong(A, C, P, invD);

	et = clock();
	//fprintf(stderr, "--- Total time at this level: %f--- \n\n", double(et - st)/CLOCKS_PER_SEC);
}
