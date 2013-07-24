% create mex files
objtocreate = {'hsc_calc_prolong', 'mxGetPropertyPtr', 'hsc_sparsify', 'ldl_solve'};
objstr = ' ';
objext = '.o';
% figure out if this is Windows
[status, result] = system('set OS');
if (~status && strcmp(result, 'Windows_NT') == 0)
    objext = 'obj';
end
for k=1:length(objtocreate)
    mexstr = sprintf('mex -largeArrayDims -c %s.cpp', objtocreate{k});
    eval(mexstr);
    objstr = sprintf('%s %s.%s', objstr, objtocreate{k}, objext);
end

mextocreate = {'hsc_precond', 'hsc_sparsify_and_compensate'};
for k=1:length(mextocreate)
    mexstr = sprintf('mex -largeArrayDims %s.cpp %s', mextocreate{k}, objstr);
    eval(mexstr);
end

