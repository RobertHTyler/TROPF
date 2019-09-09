function Lalphar = build_Lalphar(nvec, tilalpr);
% function Lalphar = build_Lalphar(nvec, tilalpr);
% 
% Builds the matrix Lalphar, which represents the operator ($ L_{\tilde{\alpha}_r} $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
% tilalpr : Attenuation coefficient(s) for rotational flow  (scalar or row vector). 
%           When tilalpr is a scalar this is a Rayleigh drag coef. 
%           When tilalpr is a row vector, it contains the coefs alpha_{r,b} 
%           for b = 0 (Rayleigh drag), b = 1 (harmonic eddy viscosity), 
%           b = 2 (biharmonic eddy viscosity), etc. When tilalpr is an N by
%           length(b) matrix each row is alpha_{r,b} at SH degree n = s,
%           s+1 ....
%           
%
% Output:
% Lalphar : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

[Nj,Mj] = size(tilalpr)     ;
N       = length(nvec)      ; % number of terms in SH expansion
 
% construct vector for velocity dissipation terms:
dissrvecs = [];
for ii = 1:Mj, 
    b  = ii-1; 
    dissrvecs(:,ii) = tilalpr(:,ii).*(-nvec.*(nvec+1)).^(b);
end
dissrvecs = sum(dissrvecs,2); dissrvecs(dissrvecs==0) = eps;


Lalphar = spdiags(   sum(dissrvecs,2)   , 0,N,N)  ;     
