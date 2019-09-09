function LAi = build_LAi(nvec,tilOm, s,tilom, Lalphad, LL);
% function LAi = build_LAi(nvec,tilOm, s,tilom, Lalphad, LL);
% 
% Builds the matrix LAi, which represents the operator $ [L_A]^{-1} $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LAi     : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LAi  = spdiags(   1./( (tilom+1i*diag(Lalphad,0)).*diag(LL) - s*tilOm )  , 0,N,N) ;  
