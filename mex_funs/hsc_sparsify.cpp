#include "mex.h"
#include<iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mxGetPropertyPtr.h"
#include "multilevel.h"

using namespace std;

/*
  Use energy-based considerations to sparsify an input matrix A.
  Assumptions are made about A being M-matrix. The resulting matrix As is triangle-free
  and a set of coarse and fine nodes are selected.
 */

mxArray *lap_to_adj(mxArray *A_IN, double *excD);
void adj_to_lap(mxArray *A_IN, mxArray *B_OUT, double *excD);
void find_nei_adj(mwIndex *A_irow, mwIndex *A_jcol, double *A_pr, int row, int * nei, int& count, int *A_irow_idx);
void sort_edges_ascend(double *W, double *Wsort, int *Widx);
void sort_edges_descend(double *W, double *Wsort, int *Widx);
void find_geom_dist_2D(int *verts, double *X, double *Y, double *dist, int *distidx);

int hsc_sparsify(mxArray *A, unsigned geom_info, double *X, double *Y, double *Z, double Ge_thresh, mxArray *C, int Level, double *tsp, double *invD, unsigned multi_triangle)
{
	mwIndex *A_irow = mxGetIr(A);
	mwIndex *A_jcol = mxGetJc(A);
	double *A_pr = (double *) mxGetPr(A);
	const mwIndex *dim = mxGetDimensions(A);
	int n = (int) dim[0];
	unsigned *C_pr = (unsigned *) mxGetPr(C);
	int divx;

	// work on adjacency matrix
	mxArray *adjA = lap_to_adj(A, invD); // convert to adjacency matrix

	A_irow = mxGetIr(adjA);
	A_jcol = mxGetJc(adjA);
	A_pr = (double *) mxGetPr(adjA);

	int *rbarray = new int[n];

	// For 2D grids - set up red-black array
	if (geom_info == 1) {
		if ((Level%2) == 0) {
			divx = (int) pow((float) 2.0, Level / 2);
			for (int i = 0; i < n; i++) {
				rbarray[i] = ((int)((X[i]/divx - 1) + (Y[i] - 1)/divx))%2;
			}
		} else {
			divx = (int) pow((float) 2.0, (Level - 1) / 2);
			for (int i = 0; i < n; i++) {
				rbarray[i] = 1 - ((int)(X[i]/divx))%2;
			}
		}
	}

	clock_t tt1, tt2;

	tt1 = clock();

	// Set up a measure of how smooth each vertex is
	double *edgeratio = new double[n], mean_edgeratio = 0.0;

	for (int i = 0; i < n; i++) {
		int ctr = 0;
		double sum = 0, min_val = 1e15, max_val = 0;
		edgeratio[i] = 0;
		for (mwIndex j = A_jcol[i]; j < A_jcol[i + 1]; j++) {
			sum += A_pr[j];
			if (A_pr[j] < min_val) { 
				min_val = A_pr[j]; 
			}
			if (A_pr[j] > max_val) { 
				max_val = A_pr[j]; 
			}
		}
		edgeratio[i] = fabs(max_val - min_val)/sum;
		mean_edgeratio = mean_edgeratio + edgeratio[i];
	}
	mean_edgeratio /= n;
	if (Ge_thresh > 1e10) {
		mean_edgeratio = Ge_thresh;
	}

	for (int v1 = 0; v1 < n; v1++) {
		// Set each point to undetermined
		C_pr[v1] = 2;
	}
  
  	/*
	 For a given triangle i, with vertices j1, j2 and j3, idx[i][k1][k2], stores the index 
	 in A_pr of A(k1, k2); where k1, k2 \in {j1, j2, j3}
	 this is to quickly update the weights during the compensation step
	*/
	int idx[1000][3][3];
	double w[3], wsort[3], dist[3];
	int verts[3], vidx[3], distidx[3], widx[3];
	int v1_nei_count = 0;
	int v1_nei[1000], v1_nei_idx[1000];

	// Set first point to F
	C_pr[0] = 0;

	for (int v1 = 0; v1 < n; v1++) {
		// Skip C points
		if (C_pr[v1] != 1) {
			// find non-zero neighbors of v1
			v1_nei_count = 0;
			find_nei_adj(A_irow, A_jcol, A_pr, v1, v1_nei, v1_nei_count, v1_nei_idx);

			// Only do for neighbors
			if (v1_nei_count > 0) {
				// cycle through neighbors of v1
				for (int v2 = 0; v2 < v1_nei_count; v2++) {
					idx[0][1][0] = v1_nei_idx[v2]; // index of A[v1_nei[v2], v1]

					// find common neighbors of v1 and v1_nei[v2]
					mwIndex v1_idx = A_jcol[v1], v2_idx = A_jcol[v1_nei[v2]];

					// index of v1 as a neighbor of v1_nei[v2]
					for (int i = v2_idx; i < (int) A_jcol[v1_nei[v2] + 1]; i++) {
						if (A_irow[i] == v1) {
							idx[0][0][1] = i; 
						}
					}

					// proceed down the neighbors of v1 and v1_nei[v2] looking for common
					// indices - these are the common neighbors; also keep checking that
					// the edge between v1 and v1_nei[v2] has not been set to 0, if so -
					// nothing to be done
					while ((v1_idx < A_jcol[v1 + 1]) && (v2_idx < A_jcol[v1_nei[v2] + 1]) && (A_pr[idx[0][0][1]] > 0)) {
						// increment the smaller of the 2 indices; if they are the same,
						// increment both
						if (A_irow[v1_idx] < A_irow[v2_idx]) {
							v1_idx++;
							continue;
						}
						if (A_irow[v2_idx] < A_irow[v1_idx]) {
							v2_idx++;
							continue;
						}

						double v1_v3_edge = A_pr[v1_idx], v2_v3_edge = A_pr[v2_idx];

						//if ((A_irow[v1_idx] == A_irow[v2_idx]) && (v1_v3_edge > 0) && (v2_v3_edge > 0)) {
						if ((v1_v3_edge > 0) && (v2_v3_edge > 0)) {
							double v1_v2_edge = A_pr[idx[0][0][1]];
							// common neighbor with non-zero connection to both v1 and
							// v1_nei[v2] found; it is A_irow[v1_idx]; 
							
							int v3 = A_irow[v1_idx];

							// triangle is (v1, v1_nei[v2], v3)
							w[0] = v2_v3_edge; w[1] = v1_v3_edge; w[2] = v1_v2_edge;
							verts[0] = v1; verts[1] = v1_nei[v2]; verts[2] = v3;

							// sort edges by ascending order of strength
							sort_edges_ascend(w, wsort, widx);

							// use geometric edge for 2D grids and 
							int vtmp0, vtmp1, vtmp2, skip_flag = 0;
							//if ((wsort[2] < Ge_thresh*wsort[0]) && X) {
							if (geom_info && edgeratio[verts[0]] <= mean_edgeratio && edgeratio[verts[1]] <= mean_edgeratio && edgeratio[verts[2]] <= mean_edgeratio && ((Level < 2) || (Ge_thresh > 1e10))) {
								// sort by  geometric distance
								find_geom_dist_2D(verts, X, Y, dist, distidx);
							
								// Edge to be cut must contain v1 as one of the vertices
								// this preserves C-C connections
								vidx[0] = distidx[2]; vidx[1] = distidx[1]; vidx[2] = distidx[0];
								vtmp0 = verts[vidx[0]]; vtmp1 = verts[vidx[1]]; vtmp2 = verts[vidx[2]];
								if (vtmp2 == v1) {
									skip_flag = 1;
								} else { 
									// in 2D grid case - set according to global R/B
									C_pr[vtmp0] = rbarray[vtmp0];
									C_pr[vtmp1] = rbarray[vtmp1];
									C_pr[vtmp2] = rbarray[vtmp2];
								}
							}
							else  {
								// order vertices by strength
								vidx[0] = widx[2]; vidx[1] = widx[1]; vidx[2] = widx[0];
								vtmp0 = verts[vidx[0]]; vtmp1 = verts[vidx[1]]; vtmp2 = verts[vidx[2]];

								// Edge to be cut must contain v1 as one of the vertices
								// this preserves C-C connections
								if (vtmp2 == v1) {
									skip_flag = 1;
								} else { 
									C_pr[vtmp0] &= 1;
									C_pr[vtmp1] &= 1;
								}
							}

							if (!skip_flag) {
								// only reach here if there is a triangle to be cut
								for (mwIndex i = A_jcol[v3]; i < A_jcol[v3 + 1]; i++) {
									if (A_irow[i] == v1) {
										idx[0][0][2] = i; // index of A[v1, v3]
									}
									else if (A_irow[i] == v1_nei[v2]) {
										idx[0][1][2] = i; // index of A[v2, v3]
									}
								}
								idx[0][2][0] = v1_idx; // index of A[v3, v1]
								idx[0][2][1] = v2_idx; // index of A[v3, v2]

								if (multi_triangle) {
									// compensate ALL triangles on this edge
									int num_common = 1; // # of common neighbors v3
									//  find the common neighbors
									for (mwIndex i = A_jcol[vtmp0]; i < A_jcol[vtmp0 + 1]; i++) {
										for (mwIndex j = A_jcol[vtmp1]; j < A_jcol[vtmp1 + 1]; j++) {
											if ((A_irow[i] == A_irow[j]) && (A_irow[i] != vtmp2) && (A_pr[i] > 0) && (A_pr[j] > 0)) {
												idx[num_common][vidx[2]][vidx[0]] = i; 
												idx[num_common][vidx[2]][vidx[1]] = j; 
												int vv = A_irow[i];
												for (mwIndex k = A_jcol[vv]; k < A_jcol[vv+1]; k++) {
													if (A_irow[k] == vtmp0) {
														idx[num_common][vidx[0]][vidx[2]] = k; 
													} else if (A_irow[k] == vtmp1) {
														idx[num_common][vidx[1]][vidx[2]] = k; 
													}
												}
												//idx[num_common][vidx[0]][vidx[1]] = idx[0][vidx[0]][vidx[1]];
												//idx[num_common][vidx[1]][vidx[0]] = idx[0][vidx[1]][vidx[0]];
												num_common++;
											}
										}
									}

									// Find conductange across each path v1-v2-v3
                           if (num_common > 10) {
                              mexErrMsgTxt("num_common > 10 - not handled! Existing..");  
                           }
									double conductance[10][2], tot_cond = 0.0;
									for (int i = 0; i < num_common; i++) {
										double w1 = A_pr[idx[i][vidx[0]][vidx[2]]], w2 = A_pr[idx[i][vidx[1]][vidx[2]]];
										// Rick's formula
										conductance[i][0] = w1; 
										conductance[i][1] = w2; 
										// Conductance-based formula
										//conductance[i][0] = (w1*w2)/(w1 + w2);
										//conductance[i][1] = (w1*w2)/(w1 + w2);
										tot_cond += conductance[i][0] + conductance[i][1];
									}

									// cut common edge 
									double wt = A_pr[idx[0][vidx[0]][vidx[1]]];
									A_pr[idx[0][vidx[0]][vidx[1]]] = 0; 
									A_pr[idx[0][vidx[1]][vidx[0]]] = 0;
									
									for (int i = 0; i < num_common; i++) {
										// now we have to cut the edge verts[vidx[0]]-verts[vidx[1]] and compensate the
										// other 2 edges

										// weight proportional to conductance
										double wt1 = 2*wt*conductance[i][0]/tot_cond; 
										double wt2 = 2*wt*conductance[i][1]/tot_cond; 
										
										A_pr[idx[i][vidx[0]][vidx[2]]] += wt1; 
										A_pr[idx[i][vidx[2]][vidx[0]]] = A_pr[idx[i][vidx[0]][vidx[2]]];
										A_pr[idx[i][vidx[1]][vidx[2]]] += wt2; 
										A_pr[idx[i][vidx[2]][vidx[1]]] = A_pr[idx[i][vidx[1]][vidx[2]]];
									} 
								} else {
									// now we have to cut the edge verts[vidx[0]]-verts[vidx[1]] and compensate the
									// other 2 edges
									double wt = A_pr[idx[0][vidx[0]][vidx[1]]];

									// cut edge 
									A_pr[idx[0][vidx[0]][vidx[1]]] = 0; 
									A_pr[idx[0][vidx[1]][vidx[0]]] = 0;
	
									A_pr[idx[0][vidx[0]][vidx[2]]] += wt; 
									A_pr[idx[0][vidx[2]][vidx[0]]] = A_pr[idx[0][vidx[0]][vidx[2]]];
									A_pr[idx[0][vidx[1]][vidx[2]]] += wt; 
									A_pr[idx[0][vidx[2]][vidx[1]]] = A_pr[idx[0][vidx[1]][vidx[2]]];
								} // multi_triangle flag
							} // skip_flag condition
						} // v1_idx == v2_idx loop 

						// if here, it means that the A_irow[v1_idx] ==
						// A_irow[v2_idx] , so increment both
						v1_idx++;
						v2_idx++;
					} // end of while loop to look for common neighbors
				} // end of v2 loop	

				// set remaining undecided neighbors to be C
				for (mwIndex i = A_jcol[v1]; i < A_jcol[v1 + 1]; i++) {
					if (C_pr[A_irow[i]] == 2) {
						if (A_pr[i] > 0) {
							C_pr[A_irow[i]] = 1;
						}
					}
				}
			} // end of v1_nei_count != 0 if statement
		} // end of C_pr[v1] == 1 if statement
	} // end of v1 loop

	int num_coarse = 0, num_fine = 0, num_undet = 0;
	for (int i = 0; i < n; i++) {
		if (C_pr[i] == 2) {
			for (mwIndex j = A_jcol[i]; j < A_jcol[i + 1]; j++) {
				if ((A_pr[j] > 0) && (C_pr[A_irow[j]] == 0)) {
					C_pr[i] = 1;
					break;
				}
			}
			// if no F neighbor, we can safely set it to F
			if (C_pr[i] == 2) {
				C_pr[i] = 0;
			}
		}
	}

	for (int i = 0; i < n; i++) {
		if (C_pr[i] == 0) {
			// Ensure that all neighbors of F are C
			for (mwIndex j = A_jcol[i]; j < A_jcol[i + 1]; j++) {
				if (A_pr[j] > 0) {
					C_pr[A_irow[j]] = 1;
				}
			}
		} else {
			num_coarse++;
		}
	}

	// finally don't have too many C's
	for (int i = 0; i < n; i++) {
		int fnei = 0;
		if (C_pr[i] == 1) {
			for (mwIndex j = A_jcol[i]; j < A_jcol[i + 1]; j++) {
				if ((A_pr[j] > 0) && (C_pr[A_irow[j]] == 0)) {
					fnei = 1;
					break;
				}
			}

			// this C has only C neighbors - we can set it to F
			if (!fnei) {
				C_pr[i] = 0;
				num_coarse--;
			}
		}
	}

	adj_to_lap(adjA, A, invD);
	
	mxDestroyArray(adjA);
	delete[] rbarray;
	delete[] edgeratio;

	tt2 = clock();
	tsp[0] = double(tt2-tt1)/CLOCKS_PER_SEC;
	fprintf(stderr, "N = %d; Coarse = %d; Fine = %d; sparsify and compensate %fs\n", n, num_coarse, n - num_coarse, tsp[0]);

	return num_coarse;
}

/*
 Return the adjacency matrix of A and also the excess diagonal
*/
mxArray *lap_to_adj(mxArray *A_IN, double *excD)
{
	mwIndex *B_irow, *B_jcol, *A_jcol, *A_irow;
    mwIndex i,j, ind;
    mwSize m,n;
    double *B, *A;
           
    m = (mwSize) mxGetM(A_IN);
    n = (mwSize) mxGetN(A_IN);

    /* Get Input */
    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A = mxGetPr(A_IN);
    
    /* Create output */
    mxArray *B_OUT = mxCreateSparse(n, n, A_jcol[n]-n, mxREAL);
    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B = mxGetPr(B_OUT);
    
    for (j = 0; j <= n; j++) {
		B_jcol[j] = A_jcol[j]-j;
	}
 
    ind = 0;  
    for (j = 0; j < n; j++) {
		double s = 0;
		for (i = A_jcol[j]; i < A_jcol[j+1]; i++) { 
			if (j!= A_irow[i]) {
               B[ind] = -A[i];
			   s -= B[ind];
               B_irow[ind] = A_irow[i];
               ind = ind+1;
			} 
			else {
				s += A[i];
			}
		}
		excD[j] = s;
	}
	return B_OUT;
}   

/* 
	Find the neighbors of row assuming an adjacency matrix is passed in.
	Also return the positions in the A_irow and A_pr arrays in A_irow_idx; and the values in vals.
*/
void find_nei_adj(mwIndex *A_irow, mwIndex *A_jcol, double *A_pr, int row, int * nei, int& count, int *A_irow_idx) 
{   
    count = 0;
	for (mwIndex i = A_jcol[row]; i < A_jcol[row+1]; i++) {
		int v = A_irow[i];
		if (A_pr[i] > 0) {
			nei[count] = v;
			A_irow_idx[count] = i;
			count++;
		}
	}
}

/*
 Convert adjacency matrix to Laplacian; assumed that space is already
 allocated for B_OUT and that the only extra element in each row of B_OUT 
 is the space for the diagonal.
*/
void adj_to_lap(mxArray *A_IN, mxArray *B_OUT, double *excD)
{
	mwIndex *B_irow, *B_jcol, *A_jcol, *A_irow;
    mwIndex i, j, k;
    double *B_pr, *A_pr;

    mwSize m = (mwSize) mxGetM(A_IN), n = (mwSize) mxGetN(A_IN);
           
    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A_pr = mxGetPr(A_IN);

    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B_pr = mxGetPr(B_OUT);
    
	// While copying over, drop zeros from the matrix
	// and entries of excD are overwritten by the inverse of diagonal elements
	int diag_idx;
	k = 0;
	B_jcol[0] = 0;
    for (j = 0; j < n; j++) {
		double s = 0;
		int cross_d = 0;
		for (i = A_jcol[j]; i < A_jcol[j + 1]; i++) { if (A_pr[i] > 0) { // non-zero element
				s += A_pr[i];
				if ((A_irow[i] > j) && (cross_d == 0)) {
					// crossed diagonal - allocate space for it to be filled
					// later
					cross_d = 1; 
					diag_idx = k;
					B_irow[k] = j;
					k++;
				}
				B_pr[k] = -A_pr[i];
				B_irow[k] = A_irow[i];
				k++;
			} 
		}
		if (cross_d == 0) {
			diag_idx = k;
			B_irow[k] = j;
			k++;
		}
		B_pr[diag_idx] = s + excD[j];
		excD[j] = 1./B_pr[diag_idx];
		B_jcol[j + 1] = k;
	}
}

/*
	Given 3 vertices in verts, and their X and Y coordinates, find their geometric distance and sort the distances.
*/
void find_geom_dist_2D(int *verts, double *X, double *Y, double *dist, int *distidx)
{
	double tmp, tmp1, dist1[3];

	tmp = X[verts[1]] - X[verts[2]];
	tmp1 = Y[verts[1]] - Y[verts[2]];
	dist1[0] = tmp*tmp + tmp1*tmp1;

	tmp = X[verts[0]] - X[verts[1]];
	tmp1 = Y[verts[0]] - Y[verts[1]];
	dist1[2] = tmp*tmp + tmp1*tmp1;

	tmp = X[verts[0]] - X[verts[2]];
	tmp1 = Y[verts[0]] - Y[verts[2]];
	dist1[1] = tmp*tmp + tmp1*tmp1;

	sort_edges_descend(dist1, dist, distidx);
}

/*
	sort 3 values in ascending order and return the indices
*/
void sort_edges_ascend(double *W, double *Wsort, int *Widx)
{
	if (W[0] <= W[1]) {
		// either W[0] or W[2] is smallet
		if (W[0] <= W[2]) {
			// W[0] <= W[2] and W[0] <= W[1]
			Wsort[0] = W[0]; Wsort[1] = W[2]; Wsort[2] = W[1];
			Widx[0] = 0; Widx[1] = 2; Widx[2] = 1;
		}
		else {
			// W[2] <= W[0] and W[2] <= W[1]
			Wsort[0] = W[2]; Wsort[1] = W[0]; Wsort[2] = W[1];
			Widx[0] = 2; Widx[1] = 0; Widx[2] = 1;
		}
	}
	else {
		// either W[1] or W[2] is smallest
		if (W[2] <= W[1]) {
			// W[2] <= W[1] and W[2] <= W[0]
			Wsort[0] = W[2]; Wsort[1] = W[1]; Wsort[2] = W[0];
			Widx[0] = 2; Widx[1] = 1; Widx[2] = 0;
		}
		else {
			// W[1] <= W[2] and W[1] <= W[0]
			Wsort[0] = W[1]; Wsort[1] = W[2]; Wsort[2] = W[0];
			Widx[0] = 1; Widx[1] = 2; Widx[2] = 0;
		}
	}
}

/*
	sort 3 values in descending order and return the indices
*/
void sort_edges_descend(double *W, double *Wsort, int *Widx)
{
	if (W[0] >= W[1]) {
		// either W[0] or W[2] is largest
		if (W[0] >= W[2]) { 
			// W[0] >= W[2] and W[0] >= W[1]
			Wsort[0] = W[0]; Wsort[1] = W[2]; Wsort[2] = W[1];
			Widx[0] = 0; Widx[1] = 2; Widx[2] = 1;
		}
		else {
			// W[2] >= W[0] and W[2] >= W[1]
			Wsort[0] = W[2]; Wsort[1] = W[0]; Wsort[2] = W[1];
			Widx[0] = 2; Widx[1] = 0; Widx[2] = 1;
		}
	}
	else {
		// either W[1] or W[2] is largest
		if (W[2] >= W[1]) {
			// W[2] >= W[1] and W[2] >= W[0]
			Wsort[0] = W[2]; Wsort[1] = W[1]; Wsort[2] = W[0];
			Widx[0] = 2; Widx[1] = 1; Widx[2] = 0;
		}
		else {
			// W[1] >= W[2] and W[1] >= W[0]
			Wsort[0] = W[1]; Wsort[1] = W[2]; Wsort[2] = W[0];
			Widx[0] = 1; Widx[1] = 2; Widx[2] = 0;
		}
	}
}

