/* ========================================================================== */
/* === CMG/Source/Solver/ldl_solve.c    ===================================== */
/* ========================================================================== */

/* 	
 *   Solve given Cholesky factorization 
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */

/* 
 * ldl_solve
 * 
 * Input:  cholesky factorization: ldl_p
 *      :  --double-- vector: b
 *
 * Output:  --double-- vector: x 
 *
 *
 */

#include "multilevel.h"
#include "stdlib.h"

void ldl_solve (ldl_p *chol, double *b, double *x)
{
    unsigned int n;
    unsigned int j;
    double *y;
    d_matrix ld, ldT;
    int i;

    n = chol->ld.n;
    ld = chol->ld;
    ldT = chol->ldT;
    y = (double *) malloc(n*sizeof(double));

    /* initialize y */
    for (i=0;i<n;i++) 
        y[i] = (double) b[chol->p[i]];
    
    for (i=0;i<n;i++) {
        for (j=ld.ia[i]; j<(ld.ia[i+1] - 1); j++) 
           y[i] = y[i]-(ld.a[j])*y[ld.ja[j]];
     }
    
    for (i=0;i<n;i++) {
       y[i]=y[i]/(ld.a[ld.ia[i+1]-1]);
    }

    for (i=(n-2);i>=0;i--) {
       for (j=(ldT.ia[i+1]-1); j>ldT.ia[i]; j--) {
          y[i] = y[i]-(ldT.a[j])*y[ldT.ja[j]];}
    }
    

    for (i=0;i<n;i++)
        x[i]= (double) y[chol->invp[i]];
    
    free(y);
}
