function Dns = DnsFrompns(pns,Kns,tilom,  LLi,LV);

%pns =  tilom^-1 * LVi * (-LL*Dns + 1i*Kns)   ;

Dns = tilom*LLi*LV*(pns) + LLi*(Kns); 