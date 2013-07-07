% create mex files
mextocreate = {'compute_cot_laplacian'};
for k=1:length(mextocreate)
    mexstr = sprintf('mex -largeArrayDims %s.cpp', mextocreate{k});
    eval(mexstr);
end

