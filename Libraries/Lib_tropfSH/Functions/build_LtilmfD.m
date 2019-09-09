function LtilmfD = build_LtilmfD(LLi,LC,LBi,LD);
% function LtilmfD = build_LtilmfD(LLi,LD,LC,LBi);
% 
% Builds the matrix Ltilp, which represents the operator $ L_{\tilde {\mathfrak D} } $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LtilmfD      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

% nvec = nVec(N,s);
 

% Composite operator:
LtilmfD  =  LLi * ( LD  -  LC * LBi * LC );

