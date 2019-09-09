function  buildParamSpace_synchrot

% % Build parameter space of solutions for tidal response of a body in a 
% synchronously rotating orbit with respect to another body. This
% limiting case can be treated generically (using the nondimensionalized
% variables) and results rescaled for specific applications.
%
% Calculations are for each of the components normalized with respect to
% |frakG| = 'Gs' for each component. Hence, the nondimensionalized
% "pressure" variable will have unit value in the equilibrium tide limit,
% meaning it becomes conformally equivalent to the forcing potential. These
% nondimensional values can be rescaled for specific applications.
%
% Several cases are calculated below (e.g. Case Rayleigh drag) in separate
% cells beginning after the cell assigning Common parameters. In each case,
% it is only the first couple lines that are edited to create a new case
% (e.g. Case Harmonic drag). 

%% % % Startup:

clear
TROPFPATH  = '~/Desktop/TROPF';     % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;   % run startup script to set paths to TROPF libraries
workdir    = [TROPFPATH,'/Studies/SynchRotation']; % working directory
resultsdir = [workdir,'/Results'];  % directory for storing *.mat file solutions
cd(workdir)

%% % % Common parameters:


% method parameters:
tilOm = 1   ; % nondim rotation rate (\Omega / \Omega_s)
N     = 500 ; % number of terms in SH expansion

% parameter space grid:
junk=[[-4:0.01:1],[1:0.05:4]]; TILCESQ=10.^junk;  junk=[-4:0.02:6]; TILT=10.^junk;
[tilcesq_TOT,tilT_TOT] = ndgrid(TILCESQ,TILT); clear junk
% if extend grid to negative tilcesq then set switch below to 1.
sw_negtilcesq = 1;


% initialize forcing coef vectors:
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
 


%% Case Rayleigh drag:
if 1==1
%
experiment_name = 'RD';
alpfacd = 1; alpfacr = 1; alpfacp = eps; sw_SG = 0;
runcalcs_synch
%
end



%% Case Newtonian radiation:
if 1==1
%
experiment_name = 'NR';
alpfacd = eps; alpfacr = eps; alpfacp = 1;  sw_SG = 0;% 
runcalcs_synch
%
end


%% Case Rayleigh drag w self gravity:
if 1==1
%
experiment_name = 'RD_SG';
alpfacd = 1; alpfacr = 1; alpfacp = eps; sw_SG = 1;
runcalcs_synch
%
end


%% Case Rayleigh drag (divergent flow only):
if 1==1
%
experiment_name = 'RDalpr0';
alpfacd = 1; alpfacr = eps; alpfacp = eps; sw_SG = 0;
runcalcs_synch
%
end

 

%% Case biharmonic dissipation:
if 1==1
%
experiment_name = 'BH';
alpfacd = [eps -1]; alpfacr = [eps -1]; alpfacp = eps; sw_SG = 0;
runcalcs_synch
%
end





%% Case Newtonian radiation w self gravity:
if 1==1
%
experiment_name = 'NR_SG';
alpfacd = eps; alpfacr = eps; alpfacp = 1;  sw_SG = 1;% 
runcalcs_synch
%
end




%% Case Rayleigh plus Newtonian radiation:
if 1==1
%
experiment_name = 'RDNR';
alpfacd = 1; alpfacr = 1; alpfacp = 1;  sw_SG = 0;% 
runcalcs_synch
%
end







%%
function runcalcs_synch
% subfunction to run calculations for tidal response due to normalized synchronous-rotation
% librations forces. Saves results to disk; no output into workspace.
%
% Eccentricity libration G20:
%
% positive tilcesq:
nF = 2; sF = 0;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G20_',experiment_name]); 
% negative tilcesq:
if sw_negtilcesq == 1
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 0;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G20N_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
end
%
%
% Eccentricity libration (westward prop component) G22W:
%
% positive tilcesq:
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G22W_',experiment_name]); 
% negative tilcesq:
if sw_negtilcesq == 1
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G22WN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
end
%
%
% Eccentricity libration (westward prop component) G22E:
%
% positive tilcesq:
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G22E_',experiment_name]); 
% negative tilcesq:
if sw_negtilcesq == 1
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G22EN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
end
%
%
% Obliquity libration (westward prop component) G21W:
%
% positive tilcesq:
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G21W_',experiment_name]); 
% negative tilcesq:
if sw_negtilcesq == 1
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G21WN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
end
%
% Obliquity libration (westward prop component) G21E:
%
% positive tilcesq:
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G21E_',experiment_name]); 
% negative tilcesq:
if sw_negtilcesq == 1
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT);
save([resultsdir,'/G21EN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
end




end % end subfunction runcalcs_synch

%%
function [workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_Tc(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,tilom,Gns,Kns,dns,ens,  tilcesq_TOT,tilT_TOT)

[Nsp,Msp] = size(tilcesq_TOT);
s = sF; 
[nvec,Ntrunc] = nVec(N,s); % vector of degrees involved

parfor iitilcesq=1:Nsp,
    percentDone = 100*iitilcesq/Nsp,
    for iitilT=1:Msp
        tilcesq = tilcesq_TOT(iitilcesq,iitilT);
        tilalpd = alpfacd/(tilT_TOT(iitilcesq,iitilT)) * abs(tilom);
        tilalpr = alpfacr/(tilT_TOT(iitilcesq,iitilT)) * abs(tilom);
        tilalpp = alpfacp/(tilT_TOT(iitilcesq,iitilT)) * abs(tilom); 
        tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
        
        if sw_SG == 1, % modify slowness vector tilnusqns for self gravity:
            rho0 = 1000; rhob = 1610; % density of ice and ocean 
            zetansvec = (3./(2*nvec+1))*(rho0/rhob);
            betansvec3 = 1 - zetansvec;
            tilnusqns  = (1./betansvec3)*tilnusqns;
        end

       
        % % Solve for Dns, then calculate Rns and pns from the Dns solution:
        [Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
        
        % % With each solution, there are some single numbers (e.g. k2, average work,
        % etc.) that usefully characterize the response. Store these in a
        % parameter-space array:
        %
        % Average several variables (densities) over globe and nondim tidal period (2pi/tilom);
        %
        % work:
        workAv_TOT(iitilcesq,iitilT)            = sum(calWns);
        % dissipation:
        dissipationAv_TOT(iitilcesq,iitilT)     = sum(calDns);
        % kinetic energy:
        kineticEnergyAv_TOT(iitilcesq,iitilT)   = sum(calEKns);
        % potential energy:
        potentialEnergyAv_TOT(iitilcesq,iitilT) = sum(calEPns);
        % Love number at the degree(s)/order of forcing:
        knFsF_TOT(iitilcesq,iitilT)  = knFsF;
        
    end % end iitilT loop
end % end iitilcesq loop
end % end function solveParamSpace_TC

end % buildParamSpace_synchrot




