#include "mex.h"
#include<list>
using namespace std;

// matrix-vector multiplication
void Ax(mxArray *A, mxArray *x, mxArray *b, list<int> idx);

// Sparse matrix declarations
// Borrowed from CMG
typedef struct smatrix {
  double *pr;
  mwIndex *irow;
  mwIndex  *jcol;
  unsigned int n;
  unsigned int issym;
} matrix;

typedef struct dmatrix {
	double    *a;
	unsigned int *ia;
	unsigned int *ja;
	unsigned int n;
	unsigned int issym;       /* true if matrix is symmetric */
} d_matrix ;


typedef struct {
	d_matrix ld;
	d_matrix ldT;
	unsigned int *p;
	unsigned int *invp;
} ldl_p ;


typedef struct shier_level {
  int cycle; // 1 for V, > 1 for W
  unsigned int *C;
  matrix P;
  matrix A;
  double *invD;
  int fine_n, coarse_n;
  matrix R; 
  matrix Rt;
  double *r, *e, *Ae;
  ldl_p chol; // Cholesky factorization 
} hier_level;

int hsc_sparsify(mxArray *A, unsigned geom_info, double *X, double *Y, double *Z, double Ge_thresh, mxArray *C, int Level, double *tsp, double *invD, unsigned multi_triangle);
void hsc_calc_prolong(mxArray *A, mxArray *C, mxArray *P, mxArray *invD);
void ldl_solve(ldl_p *chol, double *b, double *x);
