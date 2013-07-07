#include "mex.h"
#include<iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <map>
#include <list>
using namespace std;

/*
 Given a vertex list V, and a face list F, compute a cotangent Laplacian 
 L and return it. Usage: L = compute_conformal_laplacian(F, V)
*/
double myangle(double *u, double *v);

#define V prhs[0] // vertices
#define F prhs[1] // faces
#define L plhs[0] // output Laplacian

extern "C" bool mxUnshareArray(mxArray *array_ptr, bool noDeepCopy);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int nvert = 0;

	double *V_pr = (double *) mxGetData(V), *F_pr = (double *) mxGetData(F);
 	const mwSize *F_size = mxGetDimensions(F);

	for (int i = 0; i < F_size[0]*F_size[1]; i++) {
		F_pr[i] -= 1;
		if ((int) F_pr[i] > nvert) {
			nvert = (int) F_pr[i];
		}
	}
	nvert += 1;
	fprintf(stderr, "dims are %d %d; nvert %d\n", F_size[0], F_size[1], nvert);
	int nface = F_size[1];
 
 	// Currently, maximum number of non-zeros hardocoded to 10*n - must be changed later
	L = mxCreateSparse(nvert, nvert, nvert*10, mxREAL);
	double *L_pr = (double *) mxGetData(L);
	mwIndex *L_ir = mxGetIr(L), *L_jc = mxGetJc(L);

	map<int, list<int> > ring;
	// compute the faces adjacent to each vertex
	for (int i = 0; i < nface; i++) {
		for (int k = 0; k < 3; k++) {
			int vertex = F_pr[i*3 + k];
			ring[vertex].push_back(i);
		}
	}

	int nz = 0; 
	for (int i = 0; i < nvert; i++) {
		// sorted entries of row i
		// first entry is neighbor and second entry is value
		map<int, double> row;
		row.clear();
		list<int>::iterator it;

		if (i == 0) {
			it = ring[i].begin();
			while (it != ring[i].end()) {
				it++;
			}
		}
			
		ring[i].sort();
		if (i == 0) {
			it = ring[i].begin();
			while (it != ring[i].end()) {
				fprintf(stderr, "%d ", *it);
				it++;
			}
			fprintf(stderr, "\n");
		}
		
		it = ring[i].begin();
		while (it != ring[i].end()) {
			int faces[3];
			for (int k = 0; k< 3; k++) {
				faces[k] = F_pr[(*it)*3 + k];
			}
			it++;
			int vertex[2];

			if (faces[0] == i) {
				vertex[0] = faces[1]; vertex[1] = faces[2];
			} else if (faces[1] == i) {
				vertex[0] = faces[0]; vertex[1] = faces[2];
			} else if (faces[2] == i) {
				vertex[0] = faces[0]; vertex[1] = faces[1];
			} else {
				mexErrMsgTxt("Problem in face ring");
			}

			double vi[3], vj[3], vk[3];
			for (int k = 0; k < 3; k++) {
				vi[k] = V_pr[i*3 + k];		
				vj[k] = V_pr[vertex[0]*3 + k];		
				vk[k] = V_pr[vertex[1]*3 + k];		
			}

			double vki[3], vkj[3], vji[3], vjk[3];
			for (int k = 0; k < 3; k++) {
				vki[k] = vk[k] - vi[k];
				vkj[k] = vk[k] - vj[k];
				vji[k] = vj[k] - vi[k];
				vjk[k] = vj[k] - vk[k];
			}
			double alpha = myangle(vki, vkj);
			double beta = myangle(vji, vjk);

			row[vertex[0]] += 1/tan(alpha);
			row[vertex[1]] += 1/tan(beta);
		}

		// now add entries to matrix
		L_jc[i] = nz;
		map<int, double>::iterator it1 = row.begin(); 
		while (it1 != row.end()) {
			L_ir[nz] = it1->first;
			L_pr[nz] = it1->second;
			nz++;
			it1++;
		}
	}
	L_jc[nvert] = nz;


	for (int i = 0; i < F_size[0]*F_size[1]; i++) {
		F_pr[i] += 1;
	}
}

/*
 Compute angle between two 3x1 vectors
*/
double myangle(double *u, double *v)
{
	double du = 0, dv = 0, uv = 0;
	for (int k = 0; k < 3; k++) {
		du += u[k]*u[k];
		dv += v[k]*v[k];
		uv += u[k]*v[k];
	}
	du = sqrt(du);
	dv = sqrt(dv);
	if (du < 1e-16) {
		du = 1e-16;
	}
	if (dv < 1e-16) {
		dv = 1e-16;
	}
	return acos(uv/(du * dv));
}
