function Pnsx = Pnsx(n,s,x)
% Calculate Associated Legendre function Pnsx of degree=n, order=s, on x =
% cos(colatitude), where x is a scalar or a row vector (for multiple
% colatitude points)
%
% In: 
% n = degree     % scalar
% s = order      % scalar
% x = cos(colat) % may be a row vector
%
% Out:
% Pnsx = P_n^s(x) % row vector (if x is) 
%
% R. Tyler (Dec. 1, 2018)

if abs(s)>n, 
    % the forbidden case:
    Pnsx = zeros(length(x),1); % the legendre function does not allow abs(s)>n so this is generalized here 
else % the allowed cases:
    Pnsx = legendre(n,x); 
end

% extract only the associated legendre function of order/rank = s:
Pnsx=squeeze(Pnsx(s+1,:,:)); 
