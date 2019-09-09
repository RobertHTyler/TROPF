function [Dns,LtilmfD] = response_Dns(tilom,Gns,Kns,dns,ens, LVi, LL,LLi,LD,LC,LBi);
%

% build matrix representing full PDE operator for tilmfD:
LtilmfD = build_LtilmfD(N,tilOm, tilom,s,  Lalphad,Lalphar,LV);

% Solve for primary solution Dns through matrix inversion :
Dns = LtilmfD \ ( -Gns -tilom^-1*LVi*(1i*Kns) + LLi*dns + LLi*LC*LBi*ens ); 
