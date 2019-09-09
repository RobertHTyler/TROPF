function [pns,Ltilp] = response_pns(tilom,Gns,Kns,dns,ens, LVi, LL,LLi,LD,LC,LBi);
 



% % Solve for Dns, then calculate Rns and pns from the Dns solution:
Ltilp   =  LLi * ( LA - LC * LBi * LC ) * tilom * LLi * LV  + speye(N,N);
Qtilpns =  ( -Gns -tilom^-1*LVi*(1i*Kns) + LLi*dns + LLi*LC*LBi*ens ); 
pns     = (-Qtilpns)\Ltilp;

