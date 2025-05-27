function [a, r, A] = wght_spectral_fit(s,w,F) 

wF = (1./w').*F;
ws = (1./w).*s;

A = wF'*wF;
v = wF'*ws';

a = A\v;

r = s - a'*F';