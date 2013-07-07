function [Y X] = grid_2d_setup(N, M)

if (isempty(M))
  M = N;
end;

[X Y] = ind2sub([N M],[1:N*M]) ;

X = X';
Y = Y';
