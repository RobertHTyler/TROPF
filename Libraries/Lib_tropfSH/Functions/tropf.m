function [Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns)
% function [Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns)
%
% This is TROPF's primary macro function for calculating the tidal response (spherical-harmonic coefficients) 
%  as well as time/globe averages of several products.
%
% % % INPUT:
% 
% Method Parameters:
% N         : Number of terms in SH expansion
% tilOm     : Nondimensional rotation rate (usually = 1, unless e.g. rotation = 0)
%    
% Force Parameters:
% tilom:    : Temporal frequency (real scalar; positive for prograde propagation, negative for retrograde propagation)
% s:        : Order/rank of SH terms 
% Gns:      : Column vec (length N) of SH coefs for prescribed potential 
% Kns:      : Column vec (length N) of SH coefs for source/sink term in vertical structure 
% dns:      : Column vec (length N) of SH coefs for prescribed divergence of F^(p) term 
% ens:      : Column vec (length N) of SH coefs for prescribed curl of F^(p) term 
%
% Fluid Parameters:
% tilalpd   : Attenuation coefficient(s) for divergent  flow (real; scalar, or row vec (width equal number of coefs))
% tilalpr   : Attenuation coefficient(s) for rotational flow (real; scalar, or row vec (width equal number of coefs))
% tilnusqns : Squared slowness parameter (complex; scalar, or column vec (length N), if slowness varies with degree)
%
%
% % % OUTPUT:
% Dns:      : Column vec (length N) of SH coefs for the Helmholtz potential of divergent flow  
% Rns:      : Column vec (length N) of SH coefs for the Helmholtz potential of rotational flow 
% pns:      : Column vec (length N) of SH coefs for the dynamic pressure variable
% Wns:      : Column vec (length N) of globe/time averaged work rate (work at each degree; sum vec for total).
% Dns:      : Column vec (length N) of globe/time-averaged dissipation rate (dissipation at each degree; sum vec for total).
% calEKns   : Column vec (length N) of globe/time-averaged kinetic energy (kinetic energy at each degree; sum vec for total).
% calEPns   : Column vec (length N) of globe/time-averaged potential energy (potential energy at each degree; sum vec for total).
% knFsF     : Complex scalar of the admittance (Love number at the degree and order (nF,sF) of the forcing)
%
%
% % % NOTES:
%
% 1) All parameters, variables, and operators have been
% non-dimensionalized. The assumption tilOm = 1 above is equivalent to
% setting the arbitrary scaling frequency Omega_s equal the rotation rate
% Omega.
%
% 2) Abbreviations in documentation here: vec (vec); SH (spherical
% harmonic); coef (coefficient).
%
% 3) Note this is a convenient macro and is fast. But for increased speed,
% where not all these solution variables are needed, one can run the
% response*.m functions directly. For example, if one wants to calculate
% just the power---but for 10's of millions of scenarios---speed becomes
% important. One only needs Dns (or pns) to calculate power and this can be
% done with response_Dns.m or response_pns.m. That said, a profile of
% tropf.m execution shows that most of the computational time is spent on
% building the sparse matrix operators (L*) and so this macro can in fact
% save time (if multiple variables are needed) because the operators don't
% have to be rebuilt (as they would if Dns, Rns, pns... are calculated
% sequentially using the response*.m functions).
%
% 4) The variable names in the code are chosen in an attempt to reflect the
% typeset variables in the TROPF manual: "til" refers to tilde; "cal"
% refers to the caligraphic font; and "n" "s" refer to sub and superscript
% degree and order.
% 
% 10 Jan, 2019; R. Tyler %


% % Get vec of SH degrees (n = s, s+1, s+2 ... Ntrunc):
[nvec,Ntrunc] = nVec(N,s); % nvec is the vector of degrees and Ntrunc is the truncation degree.


% % Build operator matrices L*:
%
% Dissipation and slowness operators:
Lalphad   =  build_Lalphad(nvec, tilalpd);
Lalphar   =  build_Lalphar(nvec, tilalpr); 
LV        =  build_LV_fromSlowness(N,tilnusqns);
%
% Other operators:
LL        =  build_LL(nvec);
LC        =  build_LC(nvec,tilOm, s);
LD        =  build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
LVi       =  build_LVi_fromSlowness(N,tilnusqns);
LLi       =  build_LLi(nvec);
LBi       =  build_LBi(nvec,tilOm, s,tilom, Lalphar,LL);


% % Solve (one of the methods below should be uncommented):
%
% Solve for Dns, then calculate Rns and pns from the Dns solution:
LtilmfD   =  build_LtilmfD(LLi,LC,LBi,LD) ;
QtilmfD   =  build_QtilmfD(tilom,Kns,dns,ens, LVi, LLi,LC,LBi);
Dns       =  LtilmfD \ (Gns + QtilmfD) ;
Rns       =  RnsFromDns(Dns,ens, LBi,LC);
pns       =  pnsFromDns(Dns,Kns,tilom,  LVi,LL);
%
% % % Solve for pns, then calculate Dns and Rns from the pns solution:
% Ltilp     = build_Ltilp(tilom, LV, LLi,LA,LC,LBi)  ;
% Qtilp     = build_Qtilp(Kns,dns,ens,LLi,LA,LBi,LC) ;
% pns       = Ltilp \ (Gns + Qtilp)                  ; 
% Dns       = DnsFrompns(pns,Kns,tilom,  LLi,LV);
% Rns       = RnsFrompns(pns, tilom,Kns,ens, LLi,LV,LBi,LC);




% % Calculate some globe/time averaged quantities:
%
% Work rate density:
calWns  = globeTimeAverage( (-1i*Gns)             , ((LV)*(-1i*tilom*(-1i*pns))) , s ) + ...
          globeTimeAverage( ((-1i*pns)-(-1i*Gns)) , (Kns)                        , s )     ;
%
% Dissipation rate density:
calDns  = (-1/2) * globeTimeAverage( (Dns)                 , (LL*Lalphad*Dns)       , s ) ...
        + (-1/2) * globeTimeAverage( (Lalphad*Dns)         , (LL*Dns)               , s ) ...
        + (-1/2) * globeTimeAverage( (-1i*Rns)             , (LL*Lalphar*(-1i*Rns)) , s ) ...
        + (-1/2) * globeTimeAverage( (Lalphar*(-1i*Rns))   , (LL*(-1i*Rns))         , s ) ...
        + (  1 ) * globeTimeAverage( (tilom*(-1i*pns))     , (imag(LV)*(-1i*pns))   , s )   ;
%
% Kinetic Energy density:
calEKns = (-1/2) * globeTimeAverage( (Dns)     , (LL*Dns)       , s ) ... 
        + (-1/2) * globeTimeAverage( (-1i*Rns) , (LL*(-1i*Rns)) , s )   ;
%
% Potential Energy density:
calEPns = (1/2) * globeTimeAverage( (-1i*pns) , real(LV)*(-1i*pns) , s ) ;  
%
% Love number at the degree(s)/order of Gns forcing:
sF     = s;                             % order of Gns 
nF     = find(Gns) + sF - 1;            % degree(s) of non-zero Gns
knFsF  = pns((nF-sF)+1)/Gns((nF-sF)+1); % Love number at degree(s) nF

end





