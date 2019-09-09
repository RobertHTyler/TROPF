function Qtilp = build_Qtilp(Kns,dns,ens,LLi,LA,LBi,LC);


Qtilp =  - LLi*(LA - LC*LBi*LC)*LLi*(Kns) + LLi*dns + LLi*LC*LBi*ens ;
