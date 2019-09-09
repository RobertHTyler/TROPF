function LD = build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
% function LD = build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
% 
% Builds the matrix LD, which represents the operator $ L_D $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LD      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

N   = length(nvec) ; % number of terms in SH expansion


LD  = spdiags(  (  (tilom+1i*diag(Lalphad,0)).*diag(LL) - s*tilOm  ...
    + (tilom^-1)*diag(LL).*(1./diag(LV)).*diag(LL)  )          , 0,N,N) ;  

