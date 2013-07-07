function [hsc_fun hsc_hier] = hsc_hierarchy(A, opts, A_orig)

	if (strcmp(opts.legend, 'ABF') || strcmp(opts.legend, 'GMG'))
		is_abf = 1;
		if (strcmp(opts.legend, 'ABF'))
			fprintf('Switching to ABF mode!\n');
		else 
			fprintf('Switching to GMG mode!\n');
			% for GMG - set all off-diagonals to 1's
			Aorig = A;
			A = A - spdiags(diag(A), 0, size(A, 1), size(A, 2));
			[i, j, v] = find(A);
			v(:) = -1;
			A = sparse(i, j, v, size(A, 1), size(A, 2));
			sA = sum(Aorig, 2) - sum(A, 2); % make sure to add the excess diagonal back
			A = A + spdiags(sA, 0, size(A, 1), size(A, 2));
		end;
	else
		is_abf = 0;
	end;

	% threshold to switch between choosing the geometric edge and weakest edge 
	if (is_abf)
		ge_thresh = Inf;
		multi_triangle = 1;
	else
		ge_thresh = 4.0;
		multi_triangle = 0;
	end;

	t1 = clock;

	if (is_abf == 0)
		v = diag(A) - sum(A, 2);
		ratv = max(v)/min(v);
		if (ratv < 200)
			fprintf('Setting threshold to Infininity !!\n');
			ge_thresh = Inf;
			multi_triangle = 0;
		end;
	end;

	if (opts.m == 0 || opts.n == 0) 
		% 3D meshes - no grid assumed
		X = 0; 
		Y = 0; 
		geom_info = 0;
	else
		% 2D grid
		[X Y] = grid_2d_setup(opts.m, opts.n);
		geom_info = 1; 
	end;

	% check if Laplacian
	sA = sum(A, 2);
	smaxA = max(sA);
	sminA = min(sA);
	if (((abs(sminA)/abs(smaxA) ) < 1e-15) && smaxA < 1e-7) 
		fprintf('The system is numerically ill-conditioned\n');
		hsc_fun = []; 
		hsc_hier = [];
		return;
	end;

	% flag for Laplacian
	if (abs(smaxA) < 1e-4)
		lap = 1;
	else
		lap = 0;
	end;

  
	t12 = 0;
	tsp = 0;
	
	level = 1;
	hsc_hier{1}.A = A;
	hsc_hier{1}.cycle = opts.cycle;

	% currently, hardcoded to Gauss-Seidel
	hsc_hier{1}.smoother = 0; 
	hsc_hier{1}.is_laplacian = lap;
	
	hsc_hier{1}.jac_pre = opts.pre_smooth_iter;
	hsc_hier{1}.jac_post = opts.post_smooth_iter;
	hsc_hier{1}.diag_precond = opts.diag_precond;
	hsc_hier{1}.smooth_all = (1 - opts.split_spaces);

	while (level <= opts.max_levels && size(hsc_hier{level}.A, 1)>opts.direct)
		fprintf('Level %d unknowns %d nonzeros %d avg. bw %f\n', level, size(hsc_hier{level}.A,1), nnz(hsc_hier{level}.A), nnz(hsc_hier{level}.A)/ size(hsc_hier{level}.A,1));

	    n = size(hsc_hier{level}.A, 1);
		
    	% sparsify, and find interpolant and coarse variable 
		% make dummy copy so that in-place edits work in the MeX files
		As = hsc_hier{level}.A;
    	[As hsc_hier{level}.P hsc_hier{level}.C invD tsp_level] = hsc_sparsify_and_compensate(As, X, Y, [], level, ge_thresh, geom_info, multi_triangle);

		hsc_hier{level}.invD = 1./diag(hsc_hier{level}.A);

		t2 = clock;
		% Galerkin condition
    	hsc_hier{level+1}.A = (hsc_hier{level}.P'*As)*hsc_hier{level}.P;
		t12 = t12 + etime(clock, t2);
		tsp = tsp + tsp_level;

    	hsc_hier{level}.cycle = opts.cycle;
    	hsc_hier{level+1}.jac_pre = opts.pre_smooth_iter;
    	hsc_hier{level+1}.jac_post = opts.post_smooth_iter;
    	hsc_hier{level+1}.diag_precond = opts.diag_precond;
    	hsc_hier{level+1}.smooth_all = (1 - opts.split_spaces);
		hsc_hier{level+1}.is_laplacian = lap;
    
		if (geom_info) 
			coarse_idx = find(hsc_hier{level}.C == 1);
    		X = X(coarse_idx);
    		Y = Y(coarse_idx);
		end;
    	level = level+1;
  	end;

	% Coarsest level direct solver setup
	AA = hsc_hier{level}.A;

	% code borrowed from CMG
	prm = amd(AA);
	[q, d, ld] = cholGsparse(AA(prm, prm), length(AA));
	ld = ld - eye(length(d)) + diag(d);
	hsc_hier{level}.chol.ld = process_sparse_matrix(ld, 0, 2);
	hsc_hier{level}.chol.ldT = process_sparse_matrix(ld', 0, 2);
	hsc_hier{level}.chol.p = uint32(prm - 1);
	hsc_hier{level}.chol.invp = uint32(inv_permutation(prm) - 1);

	% dummy variables - required for mex functions to run correctly
	hsc_hier{level}.cycle = opts.cycle;
	hsc_hier{level}.smoother = 0;
	hsc_hier{level}.jac_pre = opts.pre_smooth_iter;
	hsc_hier{level}.jac_post = opts.post_smooth_iter;
	hsc_hier{level}.diag_precond = opts.diag_precond;
	hsc_hier{level}.smooth_all = ~opts.split_spaces;
	hsc_hier{level}.is_laplacian = lap;

   % if an A_orig is supplied use that at the finest level
   if (~isempty(A_orig))
      hsc_hier{1}.A = A_orig;
      hsc_hier{1}.invD = 1./diag(hsc_hier{1}.A);
   end;

	tt = etime(clock, t1);
	fprintf('\nHSC hierarchy setup in %f s.; Galerkin %fs; Triangle sp. %fs\n', tt, t12, tsp);
	fprintf('Size of smallest level %d\n', size(hsc_hier{level}.A,1));

	% if linear_solver_flag is set, then we want a preconditioner from only the
	% finest level - otherwise we want a preconditioner at every level.
	if (opts.linear_solver_flag)
		hsc_fun = hsc_precond_setup(hsc_hier, 1, 1.0);
	else % used for eigensolvers - so preconditioner at every level is required
		hsc_fun = hsc_precond_setup(hsc_hier, [], 1.0);
	end;
