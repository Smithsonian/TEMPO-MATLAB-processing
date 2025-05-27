function p = bilinear_warp(x,y,d)

A = zeros(3,3);
A(1,1) = length(d);
A(2,1) = sum(x);
A(3,1) = sum(y);
A(2,2) = sum(x.^2);
A(3,2) = sum(x.*y);
A(3,3) = sum(y.^2);
A(1,2) = A(2,1);
A(1,3) = A(3,1);
A(2,3) = A(3,2);

v = zeros(3,1);
v(1) = sum(d);
v(2) = sum(x.*d);
v(3) = sum(y.*d);

p = A\v;