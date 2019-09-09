function LL = build_LL(nvec);
% function LL = build_LL(nvec);
% 
% Builds the matrix LL, which represents the Laplacian operator $ L_L $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LL      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N       = length(nvec)      ; % number of terms in SH expansion
lvec    = (-nvec.*(nvec+1)) ; % vector for Laplacian coefs
lvec(lvec == 0) = eps;      ; % compensate for singular (s=0) case (otherwise Lapl inverse matrices will be singular) 

LL  = spdiags(lvec , 0,N,N) ;     

