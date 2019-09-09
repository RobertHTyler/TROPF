function LLi = build_LLi(nvec);
% function LLi = build_LLi(nvec);
% 
% Builds the matrix LLi, which represents operator $ [L_L]^{-1} $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LLi     : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N       = length(nvec)      ; % number of terms in SH expansion
lvec    = (-nvec.*(nvec+1)) ; % vector for Laplacian coefs
lvec(lvec == 0) = eps;      ; % compensate for singular (s=0) case (otherwise Lapl inverse matrices will be singular) 

LLi  = spdiags(1./lvec , 0,N,N) ;     

