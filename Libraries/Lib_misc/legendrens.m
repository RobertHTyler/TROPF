function PP = legendrens(n,s,x)
% Associated legendre function
% n = degree
% s = order
% x = cos(colat)
P=legendre(n,x);
PP=squeeze(P(s+1,:,:));
 
