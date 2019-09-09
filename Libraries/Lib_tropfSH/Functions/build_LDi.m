function LDi = build_LDi(nvec,tilOm, s,tilom, Lalphad,LV, LL);
% function LDi = build_LDi(nvec,tilOm, s,tilom, Lalphad,LV, LL);
% 
% Builds the matrix LDi, which represents the operator $ [L_D]^{-1} $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LDi     : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LDi  = spdiags(  1./(  (tilom+1i*diag(Lalphad,0)).*diag(LL) - s*tilOm  ...
    + (tilom^-1)*diag(LL).*(1./diag(LV)).*diag(LL)  )          , 0,N,N) ;  
