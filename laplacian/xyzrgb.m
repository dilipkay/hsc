function [I r g b] = xyzrgb(x,y,z)

r = 1-x ;
g = 1-y ;
b = 1-z;

I(:,:,1) = r ;
I(:,:,2) = g ;
I(:,:,3) = b ;
