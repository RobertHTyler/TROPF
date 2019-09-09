function Pnsx = Pns_xvec(n,s,x)
% Calculate Associated Legendre function Pnsx of degree=n, order=s, on x =
% cos(colatitude), where x is a scalar or a column vector (for multiple
% colatitude points)
%
% In: 
% n = degree     % scalar
% s = order      % scalar
% x = cos(colat) % scalar or column vector
%
% Out:
% Pnsx = P_n^s(x) % scalar or column vector 
%
% R. Tyler (Jan. 10, 2019)

if abs(s)>n, 
    % the forbidden case:
    Pnsx = zeros(length(x),1); % the legendre function does not allow abs(s)>n so this is generalized here 
else % the allowed cases:
    Pnsx = legendre(n,x); 
end

% extract only the associated legendre function of order/rank = s:
Pnsx = permute(squeeze(Pnsx(s+1,:,:)),[2 1]); 
