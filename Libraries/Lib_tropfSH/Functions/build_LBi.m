function LBi = build_LBi(nvec,tilOm, s,tilom, Lalphar,LL);
% function LBi = build_LBi(nvec,tilOm, s,tilom, Lalphar,LL);
% 
% Builds the inverse matrix LBi, which represents the operator ($ [L_B]^{-1} $)
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LBi     : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LBi  = spdiags(   1./((tilom+1i*diag(Lalphar,0)) .* diag(LL) - s*tilOm)  , 0,N,N) ;  
