function QtilmfD = build_QtilmfD(tilom,Kns,dns,ens, LVi, LLi,LC,LBi);


QtilmfD   =  (1/tilom)*LVi*(Kns) + LLi*dns + LLi*LC*LBi*ens ; 
