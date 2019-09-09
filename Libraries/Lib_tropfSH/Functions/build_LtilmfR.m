function LtilmfR = build_LtilmfR(N,tilOm,  Lalphad,Lalphar,LV,  tilom,s);
% function LtilmfR = build_LtilmfR(nvec,s,tilom, tilOm, LL,Lalphad,LV);
% Builds the matrix Ltilp, which represents the operator $ L_{\tilde {\mathfrak R} } $
%
% Input: 
% nvec    : Column vector of the SH degrees involved 
%           (usually nvec = nVec(n,s) for span  n = s, s+1, s+2...Ntrunc)
%
% Output:
% LtilmfR      : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018

nvec = nVec(N,s);


% Some sub operators:
LL   =  build_LL(nvec); 
LLi  =  build_LLi(nvec);
LD   =  build_LA(nvec,s,tilom, tilOm, LL,Lalphad);
LC   =  build_LC(nvec,s, tilOm);
LB   =  build_LBi(nvec,s,tilom, tilOm, LL,Lalphar);
 

LtilmfR = LLi * LC - LLi * LD * inv(LC) * LB  ; 