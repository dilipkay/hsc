function demo_3D
% demo_3D: sets up and runs a linear solver for a
% 3D mesh with cotangent Laplacian

addpath('mex_funs');
addpath('laplacian');
addpath(genpath('graph_toolbox'));

% reading is slow for large meshes
[V F] = read_off('data/lion.off');
fprintf('Mesh size %d vertices %d faces\n', size(V, 2), size(F, 2));

% add noise to the mesh 
nverts = size(V, 2);
normals = compute_normal(V, F);
sigma = 0.1; % noise level
rho = randn(1, nverts)*sigma;
V_noisy = V + (ones(3, 1)*rho).*normals;

% step size (regularization) for denoising 
% increasing this gives more denoising but
% may misshape the resulting mesh
t = 0.3;

% compute cotangent Laplacian
L = compute_cot_laplacian(V, F, 0);

% add step size
% so we are solving (tL + I) V_denoised = V_noisy
% equivalent to (L + (1/t)I) V_denoised = V_noisy/t
L = L + (1/t)*speye(size(L, 1));

% create truncated Laplacian by removing positive off-diagonals
Lt = L - spdiags(diag(L), 0, size(L, 1), size(L, 2));
[i j v] = find(Lt);
v(v > 0) = 0;
Lt = sparse(i, j, v);
% add back excess diagonal!
Lt = Lt - spdiags(sum(Lt, 2) - 1/t, 0, size(L, 1), size(L, 2)); 

% setup preconditioner on the truncated Laplacian 
hsc_fun = hsc_setup(Lt, L);

t1 = clock;
for i = 1:3
    x = pcg(L, V_noisy(i, :)'/t, 1e-6, 100, hsc_fun, []);
    V_denoised(i, :) = x';
end
fprintf('HSC solve time %g\n', etime(clock, t1));

% direct solver time taken
t1 = clock;
if (0)
   % use backslash
   for i = 1:3
      x = L\(V_noisy(i, :)'/t);
      V_denoised(i, :) = x';
   end
   fprintf('Direct solve time %g\n', etime(clock, t1));
else
   % use Cholesky - probably faster than direct solver for 
   % small meshes but quickly becomes much slower
   Chol = [];
   [Chol.L p Chol.S] = chol(L, 'lower');
   Chol.U = Chol.L';
   Chol.P = Chol.S';
   Chol.Q = Chol.S;
   fprintf('Cholesky Setup time %g\n', etime(clock, t1));
   t2 = clock;
   for i = 1:3
      x = Chol.Q * (Chol.U \ (Chol.L \ (Chol.P * (V_noisy(i, :)'/t))));
      V_denoised(i, :) = x';
   end;
   fprintf('Cholesky solve time %g\n', etime(clock, t2));
end;

% compute condition number
% slow for large systems as it requires hundreds of iterations
if (size(L, 1) < 1e5)
    fprintf('Computing Condition Number...\n');
    [c emin emax pc pemin pemax] = compute_cn(L, hsc_fun);
    fprintf('Original Condition Number %f; preconditioned system %f\n', c, pc);

    % plot
    % plotting is very slow for large meshes
    X = V(1, :);
    Y = V(2, :);
    Z = V(3, :);
    X = X - min(X); X = X./max(X);
    Y = Y - min(Y); Y = Y./max(Y);
    Z = Z - min(Z); Z = Z./max(Z);
    [I r g b] = xyzrgb(X, Y, Z);
    opts.face_vertex_color = g'/10;
    figure; plot_mesh(V, F, opts); title('Original mesh'); view([100, 0]);
    figure; plot_mesh(V_noisy, F, opts); title('Noisy Mesh'); view([100, 0]);
    figure; plot_mesh(V_denoised, F, opts); title('Denoised Mesh'); view([100, 0]);
end;
