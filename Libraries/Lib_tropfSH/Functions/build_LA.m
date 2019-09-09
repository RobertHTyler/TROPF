function LA = build_LA(nvec,tilOm, s,tilom, Lalphad, LL);
% function LA = build_LA(nvec,tilOm, s,tilom, Lalphad, LL);
% 
% Builds the matrix LA, which represents the operator ($ L_A $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LA      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LA  = spdiags(   (tilom+1i*diag(Lalphad,0)) .* diag(LL) - s*tilOm  , 0,N,N) ;  
