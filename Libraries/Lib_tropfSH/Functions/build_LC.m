function LC = build_LC(nvec,tilOm, s);
% function LC = build_LC(nvec,tilOm, s);
% 
% Builds the matrix LC, which represents the operator ($ L_C $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
% tilOm   : The nondimensional rotation rate (usually = 1)
%           (Don't confuse with tilom--the nondimensional frequency of
%           oscillation.)
% s       : Rank/order of terms in SH terms (a positive integer)
%
% Output:
% LC      : Square matrix (sparse, bidiagonal)
% 
% R. Tyler, 28 Dec., 2018

N  = length(nvec) ; % number of terms in SH expansion

LC = spdiags(   tilOm*(-nvec.*(nvec+2).*(nvec-s+1)./(2*nvec+1))   ,-1,N,N)...
   + spdiags(   tilOm*(-(nvec-1).*(nvec+1).*(nvec+s)./(2*nvec+1)) ,+1,N,N) ;
