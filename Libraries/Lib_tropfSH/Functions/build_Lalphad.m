function Lalphad = build_Lalphad(nvec, tilalpd);
% function Lalphad = build_Lalphad(nvec, tilalpd);
% 
% Builds the matrix Lalphad, which represents the operator ($ L_{\tilde{\alpha}_d} $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
% tilalpd : Attenuation coefficient(s) for divergent flow  (scalar or row vector). 
%           When tilalpd is a scalar this is the Rayleigh drag coef. 
%           When tilalpd is a row vector, it contains the coefs alpha_{d,b} 
%           for b = 0 (Rayleigh drag), b = 1 (harmonic eddy viscosity), 
%           b = 2 (biharmonic eddy viscosity), etc. 
%           
%
% Output:
% Lalphad : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018
  

[Nj,Mj] = size(tilalpd)     ;
N       = length(nvec)      ; % number of terms in SH expansion
 
% construct vector for velocity dissipation terms:
dissdvecs = [];
for ii = 1:Mj, 
    b  = ii-1; 
    dissdvecs(:,ii) = tilalpd(:,ii).*(-nvec.*(nvec+1)).^(b);
end
dissdvecs = sum(dissdvecs,2); dissdvecs(dissdvecs==0) = eps;


Lalphad = spdiags(   sum(dissdvecs,2)   , 0,N,N)  ;     
