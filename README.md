hsc
===
This package contains MATLAB and C++ Mex code to build the preconditioner
described in the following paper:

Efficient Preconditioning of Laplacian Matrices for Computer Graphics,
Dilip Krishnan, Raanan Fattal and Rick Szeliski, SIGGRAPH 2013

(c) 2013: Dilip Krishnan, Raanan Fattal and Rick Szeliski.

Please send any bug reports, questions or comments to: dilipkay@mit.edu

-------------
Installation:
-------------
1. If you have a ZIP file, unzip the package into a directory;  or 
2. You can clone from the Git repository: github.com/dilipkay/hsc
3. Run MATLAB and cd into the created directory.
4. To build various Mex executables: build_mex

------
Usage:
------
The preconditioner provided in this package works mainly with M-matrices,
which are Laplacians with non-positive off-diagonal elements. It is 
always assumed that your matrix is symmetric positive semi-definite
(this is not checked in the code for performance reasons).

Matrices which are very close to being M-matrices, with only a few 
non-negative off-diagonal elements, may be handlded by dropping the non-negative
off-diagonals and adjusting the diagonals appropriately. This can be
used for cotangent Laplacians.

If you have a Laplacian L which is an M-matrix (e.g. a graph Laplacian), 
and a right hand side b, then to use the preconditioner with Preconditioned
Conjugate Gradient, do the following:

hsc_fun = hsc_setup(L);

and then use in PCG as:

x = pcg(L, b, tolerance, max_iter, hsc_fun, []);

demo2D.m and demo3D.m provide examples of how to use the solver
to solve a 2D colorization and 3D mesh processing problem, respectively.

The core functions that set up the HSC hierarchy are: hsc_setup and hsc_hierarchy.
The parameters of the hierarchy construction are specified in hsc_setup.m. These
may be modified, for example by adding more smoothing iterations at every level
of the hierarchy. At present, only V-cycle smoothing is supported.

------------
Directories:
------------

1. data/: Contains example Laplacians and right hand sides.

2. laplacians/: Helper functions to setup cotangent Laplacians given a
OFF or PLY file.

3. graph_toolbox: Contains Gabriel Peyre's graph manipulation toolbox, with functions
for reading in OFF and  PLY files and manipulation of meshes.

-----------------
Acknowledgements:
-----------------
Yiannis Koutis for Combinatorial Multigrid. 
Gabriel Peyre for mesh manipulation toolbox.
