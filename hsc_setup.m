function hsc_fun = hsc_setup(L, L_orig, rows, cols)

% hsc_setup: Setup the HSC preconditioner for linear
% system solves. 
%
% For 2D images with known numbers of rows and columns, call as:
% hsc_fn = hsc_setup(L, rows, columns) 
%
% For 3D meshes, call as:
% hsc_fn = hsc_setup(L) 
% or as
% hsc_fn = hsc_setup(L, L_orig)
% where L_orig is the actual cotangent Laplacian and L is a "surrogate" M-matrix
% which is derived, for example, by dropping positive off-diagonals in L_orig
% and adjusting the diagonals accordingly.
%
% L is a Laplacian. Use the functions in the "laplacian" directory to set up
% Laplacians given an image or a mesh. Once the preconditioner
% hsc_fn is set up, use it in PCG with a given right-hand size b as follows:
%
% x = pcg(L, b, tol, maxit, hsc_fn);
%

if nargin == 0
	help hsc_setup
	return;
end;
addpath('mex_funs');
if (~exist('rows', 'var') || ~exist('cols', 'var'))
	m = 0;
	n = 0;
else
	m = rows;
	n = cols;
end;

if (((size(L, 1) ~= m*n) && (m > 0) && (n > 0)) || (size(L, 2) ~= size(L, 1)))
	fprintf('Some problem with rows/cols values or matrix is not square\n');
	prec_fn = [];
	return;
end;

fprintf('Setting up HSC\n');
t1  = clock;

% If 2D grid information passed in, use it
opts.m = m; 
opts.n = n;

% size of coarsest level
opts.direct = 512; 

% number of pre-smoothing iterations
% set this on for difficult problems when poor convergence is seen
opts.pre_smooth_iter = 0; 
% post-smoothing iterations
opts.post_smooth_iter = 1; 
% diagonal preconditioning of fine-level variables
opts.diag_precond = 1;

opts.legend = 'HSC';

% V-cycle (no other mode supported at present)
opts.cycle = 1;

% set to 1 to smooth all variables; 0 to only smooth C variables
opts.split_spaces = 0;

% maximum number of levels in the hierarchy
opts.max_levels = 100;

% when set to 1 - linear solver; otherwise used for eigensolver
opts.linear_solver_flag = 1;

hsc_fun = hsc_hierarchy(L, opts, L_orig);
hsc_setup_time = etime(clock, t1);

fprintf('HSC solver setup time taken %g \n', hsc_setup_time);

