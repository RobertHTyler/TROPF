function Rns = RnsFrompns(pns, tilom,Kns,ens, LLi,LV,LBi,LC);


Dns = DnsFrompns(pns,Kns,tilom,  LLi,LV);
Rns = -LBi*(LC*Dns + ens)  ;

