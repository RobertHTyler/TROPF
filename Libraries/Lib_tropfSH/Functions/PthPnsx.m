function PthPnsx = PthPnsx(n,s,x);
%function PthPnsx = PthPnsx(n,s,x);
%
% Calculate partial derivative wrt colat of Associated Legendre function P^n_s(x):
%
% Input: n (degree); s (order/rank); x=cos(colat) 
%   Notes: x can be a vector.
%   
% Output: PthPnsx is the partial derivative wrt colat of Pnsx: ($\partial_{theta} P^n_s(x)$)
%   Notes: take care to remember derivative is wrt colat and not wrt x
%  
% R. Tyler (1 Dec. 2018)


% Use of a recursion relationship for Associated Legendre functions
% provides the following:

if s-1>=0
PthPnsx    = -(1/2)*( (n+s)*(n-s+1)*Pnsx(n,s-1,x)  - Pnsx(n,s+1,x) );
elseif s-1 < 0 % else need to convert to positive argument
    MINUSm = s-1; 
    m      = -MINUSm;
PnMINUSmx  = (-1)^(m) * (1./ratiofactorials1(n,m)) *  Pnsx(n,m,x);
PthPnsx    = -(1/2)*( (n+s)*(n-s+1)*PnMINUSmx  - Pnsx(n,s+1,x) );
end


