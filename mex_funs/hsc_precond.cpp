#include "mex.h"
#include<iostream>
#include <math.h>
#include <string.h>
#include "mxGetPropertyPtr.h"
#include "multilevel.h"
using namespace std;

/*
  Given a hierarchy emg_hier, do a multigrid preconditioner on input vector ...
      x to give Ax.
 */

#define matH prhs[0] // EMG hierarchy from MATLAB
#define matx prhs[1] // input vector
#define matAx plhs[0] // output vector
#define SL prhs[2] // starting level - 1 is finest
#define ScaleFactor prhs[3] // scaling coarse level solution - for AMG it's set to 2

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

void precond(hier_level *H, double *x, double *Ax, int num_levels, int cycle, int *jac_pre, int *jac_post, int *diag_precond, int *smooth_all, int bands, int start_level, int *is_laplacian, double scale_factor);
void jacobi(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C);
void jacobi_rb(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C, int do_all, int var_order);
void mat_vec_add(matrix A, double *x, double *Ax, int m, int n);
void mat_vec(matrix A, double *x, double *Ax, int m, int n);
void mat_vec_sym(matrix A, double *x, double *Ax, int m, int n);
void residual(matrix A, double *x, double *r, double *Ax, int n);
void mat_vec_trans(matrix A, double *x, double *Ax, int m, int n);
void exact(matrix R, matrix Rt, double *b, double *x, int n);
void mat_vec_subset(matrix A, double *x, double *Ax, int n, unsigned int *C, int flag);
void gauss_seidel(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C, int var_order);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int num_levels = (int) MAX(mxGetM(matH), mxGetN(matH));
    int cycle, jac_pre[num_levels], jac_post[num_levels], diag_precond[num_levels], smooth_all[num_levels], is_laplacian[num_levels];
    hier_level *H = (hier_level *) malloc(num_levels * sizeof(hier_level));

	int start_level;
	if (prhs[2]) {
		start_level = (int) *mxGetPr(SL);
	} else {
		start_level = 1;
	}
	start_level -= 1;

	double scale_factor = (double) *mxGetPr(ScaleFactor);
		
    
    mxArray *level;

    for (int i = start_level; i < num_levels; i++) {
		level = mxGetCell(matH, i);
		if (i == start_level) {
			cycle = (int) *mxGetPr(mxGetField(level, 0, "cycle"));
		}

		mxArray *tmp;
		if (i <= (num_levels-2)) {
			H[i].C = (unsigned int *) mxGetPr(mxGetField(level, 0, "C"));
			jac_pre[i] = (int) *mxGetPr(mxGetField(level, 0, "jac_pre"));
			jac_post[i] = (int) *mxGetPr(mxGetField(level, 0, "jac_post"));
			diag_precond[i] = (int) *mxGetPr(mxGetField(level, 0, "diag_precond"));
			smooth_all[i] = (int) *mxGetPr(mxGetField(level, 0, "smooth_all"));

			tmp = mxGetField(level, 0, "invD");

			H[i].invD = (double *) mxGetPr(tmp);
			H[i].fine_n = MAX(mxGetM(tmp), mxGetN(tmp));
      
			tmp = mxGetField(level, 0, "P");
			H[i].P.pr = (double *) mxGetPr(tmp);
			H[i].P.irow = (mwIndex *) mxGetIr(tmp);
			H[i].P.jcol = (mwIndex *) mxGetJc(tmp);
			H[i].coarse_n = MIN(mxGetM(tmp), mxGetN(tmp));
		
	  		tmp = mxGetField(level, 0, "A");

			H[i].A.pr = (double *) mxGetPr(tmp);
			H[i].A.irow = (mwIndex *) mxGetIr(tmp);
			H[i].A.jcol = (mwIndex *) mxGetJc(tmp);
      	} else {
			/*tmp = mxGetField(level, 0, "cholR");
			H[i].R.pr = (double *) mxGetPr(tmp);
			H[i].R.irow = (mwIndex *) mxGetIr(tmp);
			H[i].R.jcol = (mwIndex *) mxGetJc(tmp);
			tmp = mxGetField(level, 0, "cholRt");
			H[i].Rt.pr = (double *) mxGetPr(tmp);
			H[i].Rt.irow = (mwIndex *) mxGetIr(tmp);
			H[i].Rt.jcol = (mwIndex *) mxGetJc(tmp);
			H[i].fine_n = MAX(mxGetM(tmp), mxGetN(tmp));*/
			tmp = mxGetField(level, 0, "A");
			H[i].fine_n = MAX(mxGetM(tmp), mxGetN(tmp));
			H[i].chol.ld.a = (double *) mxGetPr(mxGetField(mxGetField(mxGetField(level , 0, "chol"), 0, "ld"), 0, "a")); 
			H[i].chol.ld.ia = (unsigned int *) mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ld"), 0, "ia"));        
			H[i].chol.ld.ja = (unsigned int *) mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ld"), 0, "ja"));       
			H[i].chol.ld.n = (unsigned int ) *(mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ld"), 0, "n")));
			H[i].chol.ldT.a = (double *) mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ldT"), 0, "a"));
			H[i].chol.ldT.ia = (unsigned int *) mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ldT"), 0, "ia")); 
			H[i].chol.ldT.ja = (unsigned int *) mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ldT"), 0, "ja"));
			H[i].chol.ldT.n = (unsigned int ) *(mxGetPr(mxGetField(mxGetField(mxGetField(level, 0, "chol"), 0, "ldT"), 0, "n"))); 
			H[i].chol.p   = (unsigned int *) mxGetPr(mxGetField(mxGetField(level, 0, "chol"), 0, "p")); 
			H[i].chol.invp   = (unsigned int *) mxGetPr(mxGetField(mxGetField(level, 0, "chol"), 0, "invp"));

			// debug
			tmp = mxGetField(level, 0, "A");
			H[i].A.pr = (double *) mxGetPr(tmp);
			H[i].A.irow = (mwIndex *) mxGetIr(tmp);
			H[i].A.jcol = (mwIndex *) mxGetJc(tmp);
      	}
		tmp = mxGetField(level, 0, "r");
		if (i > start_level) {
			H[i].r = (double *) mxGetPr(tmp);
			tmp = mxGetField(level, 0, "e");
			H[i].e = (double *) mxGetPr(tmp);
		}
		tmp = mxGetField(level, 0, "Ae");
		H[i].Ae = (double *) mxGetPr(tmp);
		is_laplacian[i] = (int) *mxGetPr(mxGetField(level, 0, "is_laplacian"));
    } 

    double *x_pr = (double *) mxGetPr(matx);
    matAx = mxCreateNumericArray(mxGetNumberOfDimensions(matx), mxGetDimensions(matx), mxDOUBLE_CLASS, mxREAL);

    double *Ax_pr = mxGetPr(matAx);

    precond(H, x_pr, Ax_pr, num_levels, cycle, jac_pre, jac_post, diag_precond, smooth_all, mxGetN(matx), start_level, is_laplacian, scale_factor);
}

/* 
 Do the actual preconditioning - by sweeping up and down the hierarchy. 
 TODO: add Gauss-Seidel smoothing and for multiple pre- and post-smoothing 
 steps and different omega.
*/

void precond(hier_level *H, double *x, double *Ax, int num_levels, int cycle, int *jac_pre, int *jac_post, int *diag_precond, int *smooth_all, int bands, int start_level, int *is_laplacian, double scale_factor)
{
	double omega = 0.8;
	double *Ae[num_levels];

	for (int b = 0; b < bands; b++) {
		H[start_level].r = &(x[b*H[start_level].fine_n]);
		H[start_level].e = &(Ax[b*H[start_level].fine_n]);

		for(int iter = 0 ; iter < cycle; iter++) {
			// downward sweep - finest to coarsest
			for(int i = start_level; i < num_levels - 1; i++) {
				int n = H[i].fine_n;

				// pre-smoothing
				if (jac_pre[i]) {
					if (iter == 0) {
						memset(H[i].e, 0, H[i].fine_n*sizeof(double));
					}
					//jacobi_rb(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_pre[i], H[i].C, smooth_all[i], 0);
					//jacobi(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_pre[i], H[i].C);
					gauss_seidel(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_pre[i], H[i].C, 0);
					Ae[i] = H[i].Ae;
					// compute residual: Ae = r - A*e
					residual(H[i].A, H[i].e, H[i].r, H[i].Ae, n);
				} else {
					if (iter == 0) {
						// no pre-smoothing: so residual = r
						//memcpy(H[i].Ae, H[i].r, H[i].fine_n*sizeof(double));
						Ae[i] = H[i].r;
					} else {
						residual(H[i].A, H[i].e, H[i].r, H[i].Ae, n);
						Ae[i] = H[i].Ae;
					}
				}
				
				// restrict
				mat_vec_trans(H[i].P, Ae[i], H[i + 1].r, H[i].fine_n, H[i].coarse_n);

				if ((iter>0) || ((jac_pre[i]) && (b > 0))) {
					memset(H[i + 1].e, 0, H[i+1].fine_n*sizeof(double));
				}
			}
		
			// exact solve at coarsest level
			// e = A^{-1} r
			// Ensure coarsest level inversion is stable by removing mean
			if (is_laplacian[start_level]) {
				double mean_r = 0.0;
				for (int i = 0; i < H[num_levels - 1].fine_n; i++) {
					mean_r += H[num_levels - 1].r[i];
				}
				mean_r /= (double) H[num_levels - 1].fine_n;
				for (int i = 0; i < H[num_levels - 1].fine_n; i++) {
					H[num_levels - 1].r[i] -= mean_r;
				}
			}
			//exact(H[num_levels-1].R, H[num_levels-1].Rt, H[num_levels-1].r, H[num_levels-1].e, H[num_levels-1].fine_n);
			ldl_solve(&(H[num_levels - 1].chol), H[num_levels - 1].r, H[num_levels - 1].e);


			if (0) {
				double acc = 0, tmpv[H[num_levels-1].fine_n], tmpv1[H[num_levels-1].fine_n];
				mat_vec(H[num_levels-1].A, H[num_levels-1].e, tmpv1, H[num_levels-1].fine_n, H[num_levels-1].fine_n);
				for (int i=0; i<H[num_levels-1].fine_n; i++) {
					acc += fabs(tmpv1[i] - H[num_levels-1].r[i]);
				}
				fprintf(stderr, "exact solve acc. %f\n", acc);
			}
			

			for(int i = num_levels-2; i >= start_level; i--) {
				// prolong error and add to that of previous level if there is
				// pre-smoothing, else just interpolate error to finer level

				if (scale_factor > 1.0) {
					for (int j = 0; j < H[i + 1].fine_n; j++) {
						H[i + 1].e[j] *= scale_factor;
					}
				}

				if (jac_pre[i]) {
					mat_vec_add(H[i].P, H[i+1].e, H[i].e, H[i].fine_n, H[i].coarse_n);
				} else {
					mat_vec(H[i].P, H[i+1].e, H[i].e, H[i].fine_n, H[i].coarse_n);
				}

				if (diag_precond[i]) {
					// do diagonal preconditioning of fine variables
					// on residual before coarse-grid correction
					for (int j=0; j<H[i].fine_n; j++) {
						if (H[i].C[j] == 0) {
							H[i].e[j] += Ae[i][j]*H[i].invD[j];
						}
					}
				}

				// post-smoothing
				if (jac_post[i]) {
					//jacobi_rb(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_post[i], H[i].C, smooth_all[i], 1);
					//jacobi(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_post[i], H[i].C);
					gauss_seidel(H[i].A, H[i].invD, H[i].r, H[i].e, omega, H[i].fine_n, jac_post[i], H[i].C, 1);
				}
			}

			/*// Clamp one variable for Laplacians
			if (is_laplacian[start_level]) {
				double mean_e = 0.0;
				for (int i = 0; i < H[start_level].fine_n; i++) {
					mean_e += H[start_level].e[i];
				}
				mean_e /= (double) H[start_level].fine_n;
				for (int i = 0; i < H[start_level].fine_n; i++) {
					H[start_level].e[i] -= mean_e;
				}
			}*/
		}
	}
}

/*
  Jacobi smoothing - currently only a single iteration. 
  Note that function assumes that x is the initial iterate
  as well as output: so x must not contain garbage.
 */
void jacobi(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C) 
{
  double *Ax = (double *) malloc(n*sizeof(double));

  for (int j=0; j<iter; j++) {
    mat_vec(A, x, Ax, n, n);

    for (int i=0; i<n; i++) 
      x[i] += w*(rhs[i] - Ax[i])*invD[i];
  }

  free(Ax);
}

/*
  Jacobi-Fine-Coarse smoothing - currently only a single iteration. 
  Note that function assumes that x is the initial iterate
  as well as output: so x must not contain garbage. 
 */
void jacobi_rb(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C, int do_all, int var_order) 
{
	// to ensure positive definiteness of the preconditioner, we have to reverse
	// variable ordering in the pre-smoothing and post-smoothing cases
	if (var_order) {
		for (int j=0; j<iter; j++) {
			if (do_all) {
				// First one set of points
				for (int i = 0; i < n; i++) {
					if (C[i] == 0) {
						double tmp = 0;
						for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
							int row = A.irow[k];
							tmp += A.pr[k]*x[row];
    					}
						x[i] += (rhs[i] - tmp)*invD[i];
					}
  				}
			}

			// then the other set 
			for (int i = 0; i < n; i++) {
				if (C[i] == 1) {
					double tmp = 0;
					for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
						int row = A.irow[k];
						tmp += A.pr[k]*x[row];
					}
					x[i] += (rhs[i] - tmp)*invD[i];
				}
			}
		}
	} else {
		for (int j=0; j<iter; j++) {
			for (int i = n - 1; i >= 0; i--) {
				if (C[i] == 1) {
					double tmp = 0;
					for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
						int row = A.irow[k];
						tmp += A.pr[k]*x[row];
					}
					x[i] += (rhs[i] - tmp)*invD[i];
				}
			}
			if (do_all) {
				for (int i = n - 1 ; i >= 0; i--) {
					if (C[i] == 0) {
						double tmp = 0;
						for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
							int row = A.irow[k];
							tmp += A.pr[k]*x[row];
    					}
						x[i] += (rhs[i] - tmp)*invD[i];
					}
  				}
			}
		}
	}
}

/*
  Gauss-Seidel smoothing - currently only a single iteration. 
  Note that function assumes that x is the initial iterate
  as well as output: so x must not contain garbage. 
 */
void gauss_seidel(matrix A, double *invD, double *rhs, double *x, double w, int n, int iter, unsigned int *C, int var_order) 
{
	// to ensure positive definiteness of the preconditioner, we have to reverse
	// variable ordering in the pre-smoothing and post-smoothing cases
	if (var_order) {
		for (int j=0; j<iter; j++) {
			// First one set of points
			for (int i = 0; i < n; i++) {
				double tmp = 0;
				for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
					int row = A.irow[k];
					tmp += A.pr[k]*x[row];
    			}
				x[i] += (rhs[i] - tmp)*invD[i];
			}
		}
	} else {
		for (int j=0; j<iter; j++) {
			for (int i = n - 1; i >= 0; i--) {
				double tmp = 0;
				for (int k = A.jcol[i]; k < A.jcol[i + 1]; k++) {
					int row = A.irow[k];
					tmp += A.pr[k]*x[row];
				}
				x[i] += (rhs[i] - tmp)*invD[i];
			}
		}
	}
}


/*
  Matrix vector multiply with overwriting of output vector
 */
void mat_vec(matrix A, double *x, double *Ax, int m, int n)
{
  for (int i=0; i<m; i++) 
    Ax[i] = 0;

  for (int i=0; i<n; i++) {
    for (int j=A.jcol[i]; j<A.jcol[i+1]; j++) {
      int row = A.irow[j];
      Ax[row] += A.pr[j]*x[i];
    }
  }
}

/*
  Faster mat-vec taking into account symmetry of A
 */
void mat_vec_sym(matrix A, double *x, double *Ax, int m, int n)
{
	/*
		Ax_i = \sum_j A_ij * x_j = \sum_j A_ji  * x_j
	*/
	for (int i = 0; i < n; i++) {
		Ax[i] = 0;	
		for (int j = A.jcol[i]; j < A.jcol[i+1]; j++) {
			int row = A.irow[j];
			Ax[i] += A.pr[j]*x[row];
		}
	}
}

/*
  Compute residual Ax = r - A*x; assumes A is symmetric
 */
void residual(matrix A, double *x, double *r, double *Ax, int n)
{
	/*
		x_i = \sum_j A_ij * x_j = \sum_j A_ji  * x_j
	*/
	for (int i = 0; i < n; i++) {
		Ax[i] = 0;	
		for (int j = A.jcol[i]; j < A.jcol[i+1]; j++) {
			int row = A.irow[j];
			Ax[i] += A.pr[j]*x[row];
		}
		Ax[i] = r[i] - Ax[i];
	}
}

/*
  Matrix vector multiply with overwriting of output vector; only do for those points
  where C[i] == flag. Assumed that A is square and symmetric
 */
void mat_vec_subset(matrix A, double *x, double *Ax, int n, unsigned int *C, int flag)
{
	for (int i = 0; i < n; i++) {
		if (C[i] != flag) {
			continue;
		}
		Ax[i] = 0;
		for (int j = A.jcol[i]; j < A.jcol[i + 1]; j++) {
			int row = A.irow[j];
			Ax[i] += A.pr[j]*x[row];
    }
  }
}


/*
  Matrix vector multiply with NO overwriting of output vector
*/

void mat_vec_add(matrix A, double *x, double *Ax, int m, int n)
{
  for (int i=0; i<n; i++) {
    for (int j=A.jcol[i]; j<A.jcol[i+1]; j++) {
      int row = A.irow[j];
      Ax[row] += A.pr[j]*x[i];
    }
  }
}

/*
  Matrix vector transpose-multiply with overwriting of output vector.
  Here A is mxn, but the output is A'x.
 */
void mat_vec_trans(matrix A, double *x, double *Atx, int m, int n)
{
  /*for (int i=0; i<n; i++) 
    Atx[i] = 0;*/

  for (int i=0; i<n; i++) {
    double val = 0;
    for (int j=A.jcol[i]; j<A.jcol[i+1]; j++) {
      val += A.pr[j]*x[A.irow[j]];
    }
    Atx[i] = val;
  }
}

/*
  Do exact solve of Ax = b using the upper triangular 
  matrix R i.e. solve R'Rx = b.
 */

void exact(matrix R, matrix Rt, double *b, double *x, int n)
{
  double *z = (double *) malloc(n*sizeof(double));

  // Solve R'z = b
  for (int i=0; i<n; i++) {
    double val = 0;
    // compute sum_{j \ne i} R'(i,j)z_j == sum_{j \ne i} R(j,i)z_j
    // in  the loop assuming that R(i,i) is the last element 
    // of the ith column of R
    for (int j=R.jcol[i]; j<R.jcol[i+1]-1; j++) {
      val += R.pr[j]*z[R.irow[j]];
    }
    z[i] = (b[i]-val)/R.pr[R.jcol[i+1]-1];
  }

  // Then Rx = z
  for (int i=n-1; i>=0; i--) {
    double val = 0;
    // compute sum_{j \ne i} R(i,j)z_j == sum_{j \ne i} Rt(j,i)z_j
    // in the loop, assuming that Rt(i,i) is the first element
    // of the i'th column of R'
    for (int j=Rt.jcol[i]+1; j<Rt.jcol[i+1]; j++) {
      val += Rt.pr[j]*x[Rt.irow[j]];
    }
    x[i] = (z[i]-val)/Rt.pr[Rt.jcol[i]];
  }
  free(z);
}
