



%% Housekeeping:
clear
TROPFPATH  = '~/Desktop/TROPF';     % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;   % run startup script to set paths to TROPF libraries
workdir    = [TROPFPATH,'/Validations'];
cd(workdir)


%% % Compare results from tropf and tropf_sas to show they give the same results (to within ~1e-18)

% Example input parameters:


% method parameters:
tilOm = 1   ; % nondim rotation rate (\Omega / \Omega_s)
N     = 500 ; % number of terms in SH expansion

% Eccentricity libration G22:
nF = 2; sF = 2;       % degree and order of forcing potential 
s = sF;               % all terms in expansion have s = sF
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = 1i*( (-1/2) / PnFsF_amp ); 
%alpfacd = 1; alpfacr = 1; alpfacp = 0;

%tilalpd   = [rand(1),rand(1),rand(1)]; tilalpr = [rand(1),rand(1)]; tilalpp = rand(1);
tilalpd = 1; tilalpr = 1; tilalpp = 0; 
tilcesq   = rand(1);
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;

 
% Solve using several methods:

% Using tropf_sas:
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns, M1,M2,F1,F2,M1e,M2e] = tropf_sas(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns); Dns_N = Dns; Rns_N = Rns; pns_N = pns; 
calWns_N = calWns; calDns_N = calDns; calEKns_N = calEKns; calEPns_N = calEPns;  

% Using tropf (this is the primary method):
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);


mabs(Dns_N,Dns)
mabs(Rns_N,Rns)
mabs(pns_N,pns)
mabs(calWns_N,calWns)
mabs(calDns_N,calDns)
mabs(calEKns_N,calEKns)
mabs(calEPns_N,calEPns)


 