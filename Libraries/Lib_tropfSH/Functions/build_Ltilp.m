function Ltilp = build_Ltilp(tilom, LV, LLi,LA,LC,LBi);
% function Ltilp = build_Ltilp(tilom, LV, LLi,LA,LC,LBi);
%
% Builds the matrix Ltilp, which represents the operator $ L_{\tilde p } $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% Ltilp   : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018
 
N = length(LLi); 

% Composite operator:
Ltilp  =  LLi * ( LA - LC * LBi * LC ) * tilom * LLi * LV  + speye(N,N);
