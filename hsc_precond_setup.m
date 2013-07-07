% Given the hierarchy hsc_hier
% set up a function hsc_fun which takes in a vector x and
% optionally a start level and preconditions the vector x from
% the start level to coarsest level
%
function hsc_fun = hsc_precond_setup(hsc_hier, start_level, scale_factor)
	for i = 1:length(hsc_hier)
		hsc_hier{i}.e = zeros(size(hsc_hier{i}.A, 1), 1);	
		hsc_hier{i}.r = zeros(size(hsc_hier{i}.A, 1), 1);	
		hsc_hier{i}.Ae = zeros(size(hsc_hier{i}.A, 1), 1);	
	end;

	if (start_level == 1)
		hsc_fun = @(x) hsc_precond(hsc_hier, x, 1, scale_factor);
	else
		hsc_fun = @(x, sl) hsc_precond(hsc_hier, x, sl, scale_factor);
	end;
