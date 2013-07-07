% create mex files
mextocreate = {'compute_conformal_laplacian'}
for k=1:length(mextocreate)
    mexstr = sprintf('mex -largeArrayDims %s.cpp', mextocreate{k});
    eval(mexstr);
end

