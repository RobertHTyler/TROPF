function sPnsxOVsinth = sPnsxOVsinth(n,s,x);
%function sPnsxOVsinth = sPnsxOVsinth(n,s,x);
%
% Calculate ratio s*Pnsx/sin(colat), where Pnsx is associated legendre function (P^n_s(x)):
%
% Input: n (degree); s (order/rank); x=cos(colat) 
%   Notes: x can be a vector.
%   
% Output: sPnsOVsinth ($\partial_{theta} P^n_s$)

% Simply division gives the following:
sinth = sqrt(1-x.^2);
sPnsxOVsinth = (s*Pnsx(n,s,x)) ./ sinth ;

