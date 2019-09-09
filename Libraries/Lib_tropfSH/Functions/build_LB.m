function LB = build_LB(nvec,tilOm, s,tilom, Lalphar, LL);
% function LB = build_LB(nvec,tilOm, s,tilom, Lalphar, LL);
% 
% Builds the matrix LB, which represents the operator ($ L_B $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LB      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LB  = spdiags(   (tilom+1i*diag(Lalphar,0)) .* diag(LL) - s*tilOm  , 0,N,N) ;  
