  function [Dns,Rns,pns, calWns,calDns,calEKns,calEPns, M1,M2,F1,F2,M1e,M2e] = tropf_sas(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns)
%function [Dns,Rns,pns, calWns,calDns,calEKns,calEPns,tiltauns, M1,M2,F1,F2,M1e,M2e] = tropf_sas(N,tilOm, tilalpd,tilalpr,tilnusqns, tilom,s,Gns,Kns,dns,ens)
%
% This is TROPF's alternative macro function for calculating the tidal
% response (spherical-harmonic coefficients) as well as time/globe averages
% of several products. This routine is formulated through solving the
% symmetric (about equator) and antisymmetric systems for the two velocity
% potentials (the "sas" in the function name stands for
% "Symmetric-AntiSymmetric"). This routine was developed and used by R.
% Tyler but has been largely obsoleted by tropf.m (which uses a different
% formulation). tropf_sas.m is used in the TROPF validations to show that
% tropf_sas.m and tropf.m produce the same results to high precision.
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
% Dns:      vec of SH coefs for the Helmholtz potential of divergent flow (i.e. \mathfrak{D_n_s})
% Rns:      vec of SH coefs for the Helmholtz potential of divergent flow (i.e. \mathfrak{R_n_s})
% pns:      : Column vec (length N) of SH coefs for the dynamic pressure variable
% Wns:      : Column vec (length N) of globe/time averaged work rate (work at each degree; sum vec for total).
% Dns:      : Column vec (length N) of globe/time-averaged dissipation rate (dissipation at each degree; sum vec for total).
% calEKns   : Column vec (length N) of globe/time-averaged kinetic energy (kinetic energy at each degree; sum vec for total).
% calEPns   : Column vec (length N) of globe/time-averaged potential energy (potential energy at each degree; sum vec for total).
% knFsF     : Complex scalar of the admittance (Love number at the degree and order (nF,sF) of the forcing)
%
% M1, M1e:  : coef matrix for symmetric system
% M2, M2e:  : coef matrix for antisymmetric system
% F1, F2:   : Column vecs (length N) of SH coefs for forcing tidal potential (symm and antisymm)
%
%
% % % NOTES:
% 1) All parameters, variables, and operators have been non-dimensionalized.
% 2) Radial symmetry in the material parameters is assumed.
% 3) Abbreviations in documentation here: vec (vector); SH (spherical harmonic); coef (coefficient)
%
% Revision history:
% 23, July, 2018; R. Tyler
%
% 
% %
%%
 
% % Get vec of degrees (n = s, s+1, s+2 ... Ntrunc):
[nvec,Ntrunc] = nVec(N,s);


% % Build operator matrices:
%
% Dissipation and slowness operators:
Lalphad   =  build_Lalphad(nvec, tilalpd);
Lalphar   =  build_Lalphar(nvec, tilalpr);
LV        =  build_LV_fromSlowness(N,tilnusqns);
%
% Other operators:
LL        =  build_LL(nvec);
LA        =  build_LA(nvec,tilOm, s,tilom, Lalphad, LL);
LB        =  build_LB(nvec,tilOm, s,tilom, Lalphar, LL);
LC        =  build_LC(nvec,tilOm, s);
LD        =  build_LD(nvec,tilOm, s,tilom, Lalphad,LV, LL);
LVi       =  build_LVi_fromSlowness(N,tilnusqns);
 


% 
M1 = sparse(N,N); 
M2 = sparse(N,N);
junk =  LD + LC;
M1(1:2:end,:) =  junk(1:2:end,:);
M2(2:2:end,:) =  junk(2:2:end,:);
junk =  LB + LC; 
M2(1:2:end,:) =  junk(1:2:end,:);
M1(2:2:end,:) =  junk(2:2:end,:);




% 
LLi       =  build_LLi(nvec);
M1e = sparse(N,N); 
M2e = sparse(N,N);
junk =  tilom*LLi*LLi * (LA + LC);
M1e(1:2:end,:) =  junk(1:2:end,:);
M2e(2:2:end,:) =  junk(2:2:end,:);
junk =   tilom*LLi*LLi * (LB + LC); 
M2e(1:2:end,:) =  junk(1:2:end,:);
M1e(2:2:end,:) =  junk(2:2:end,:);


%% %%% Now solve for SH-coef vectors Dns, Rns:

[tilcesq,tilalpp] = cNalppFromSlow(tilom,tilnusqns)
 


% Build force coef vectors F1 (for equatorially symmetric system), F2 (asymmetric system)
F1 = zeros(N,1); F2 = zeros(N,1);
Qs1ns = (1/tilom)*LL*LVi*Kns + dns;
Qs2ns = -ens; 
%
junk = LL*Gns;
F1(1:2:N) = junk(1:2:N) + Qs1ns(1:2:N);
F1(2:2:N) =               Qs2ns(2:2:N); 
F2(1:2:N) =               Qs2ns(1:2:N);
F2(2:2:N) = junk(2:2:N) + Qs1ns(2:2:N);
 

% Solution vectors:
SOL1 = M1 \ (F1);
SOL2 = M2 \ (F2);



% reshuffle to get coef vectors A and B:
Dns=zeros(N,1); Rns=zeros(N,1);
Dns(1:2:N) = SOL1(1:2:N,1);
Dns(2:2:N) = SOL2(2:2:N,1);
Rns(2:2:N) = SOL1(2:2:N,1);
Rns(1:2:N) = SOL2(1:2:N,1);


if 1==2
for ii=1:2:N
Dns(ii,1)   = SOL1(ii,1);
Rns(ii+1,1) = SOL1(ii+1,1);
end
for ii=1:2:N
Dns(ii+1,1) = SOL2(ii+1,1);
Rns(ii,1)   = SOL2(ii,1);
end
end


pns   = pnsFromDns(Dns,Kns,tilom,  LVi,LL);

 

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






