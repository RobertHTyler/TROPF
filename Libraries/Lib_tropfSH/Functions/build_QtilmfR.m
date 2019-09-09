function QtilmfR = build_QtilmfR(tilom,Kns,dns,ens, LVi, LLi,LD,LC);


QtilmfR   =  (1/tilom)*LVi*(Kns) + LLi*dns + LLi*LD*inv(LC)*ens ; 
