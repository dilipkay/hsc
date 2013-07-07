function demo_2D
% demo_2D: sets up and runs a linear solver for a
% 2D colorization problem

addpath('mex_funs');
load data/colorization.mat;

hsc_fun = hsc_setup(L, L, rows, cols);

t1 = clock;
x = pcg(L, b(:, 1), 1e-6, 100, hsc_fun, []);
fprintf('Solve time %g\n', etime(clock, t1));

% condition number computation: very slow for large meshes
fprintf('Computing Condition Number...\n');
[c emin emax pc pemin pemax] = compute_cn(L, hsc_fun);
fprintf('Original Condition Number %f; preconditioned system %f\n', c, pc);




