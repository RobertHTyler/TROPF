function  buildParamSpace_arbrot

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
workdir    = [TROPFPATH,'/Studies/NonSynchRotation/ArbRotation']; % working directory
resultsdir = [workdir,'/Results'];  % directory for storing *.mat file solutions
cd(workdir)

%% % % Common parameters:


% method parameters:
tilOm = 1   ; % nondim rotation rate (\Omega / \Omega_s)
N     = 500 ; % number of terms in SH expansion

% parameter space grid:
junk    = [[-4:0.01:1],[1+0.05:0.05:4]]; 
TILCESQ = 10.^junk;      % samples of tilcesq parameter
TILOM   = [-1:.001:1]*2; % samples of tilom parameter 
[tilcesq_TOT,tilom_TOT] = ndgrid(TILCESQ,TILOM); % grid of (tilcesq, tilom) parameters
clear junk


% define tilT:
tilT = 1e2;
 
% initialize forcing coef vectors:
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
 


%% Case Rayleigh drag:
if 1==1
%
experiment_name = 'RD';
alpfacd = 1; alpfacr = 1; alpfacp = eps; sw_SG = 0;
runcalcs_arbrot
%
end



%% Case Newtonian radiation:
if 1==1
%
experiment_name = 'NR';
alpfacd = eps; alpfacr = eps; alpfacp = 1;  sw_SG = 0;% 
runcalcs_arbrot
%
end



%% Case Rayleigh drag w self gravity:
if 1==2
%
experiment_name = 'RD_SG';
alpfacd = 1; alpfacr = 1; alpfacp = eps; sw_SG = 1;
runcalcs_arbrot
%
end


%% Case Rayleigh drag (divergent flow only):
if 1==2
%
experiment_name = 'RDalpr0';
alpfacd = 1; alpfacr = eps; alpfacp = eps; sw_SG = 0;
runcalcs_arbrot
%
end

 

%% Case biharmonic dissipation:
if 1==2
%
experiment_name = 'BH';
alpfacd = [eps -1]; alpfacr = [eps -1]; alpfacp = eps; sw_SG = 0;
runcalcs_arbrot
%
end




%% Case Newtonian radiation w self gravity:
if 1==2
%
experiment_name = 'NR_SG';
alpfacd = eps; alpfacr = eps; alpfacp = 1;  sw_SG = 1;% 
runcalcs_arbrot
%
end




%% Case Rayleigh plus Newtonian radiation:
if 1==2
%
experiment_name = 'RDNR';
alpfacd = 1; alpfacr = 1; alpfacp = 1;  sw_SG = 0;% 
runcalcs_arbrot
%
end







%%
function runcalcs_arbrot
% subfunction to run calculations for tidal response due to arbitrary
% non-synch rotation (only relative rotation term calculated here; no
% eccentricity/obliquity components)
% 
%
% Sectoral tidal force potential (G22W):
%
% positive tilcesq:
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
%tilom      = -1;      % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_cesqom(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,Gns,Kns,dns,ens,  tilT,tilcesq_TOT,tilom_TOT);
save([resultsdir,'/G22W_',experiment_name]); 
% negative tilcesq:
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 2;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
%tilom      = -1;      % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_cesqom(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,Gns,Kns,dns,ens,  tilT,tilcesq_TOT,tilom_TOT);
save([resultsdir,'/G22WN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
%
%
% Tesseral tidal force potential (G21W):
%
% positive tilcesq:
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;       % amplitude of associated Legendre function of degree nF, order sF;
%tilom      = -1/2;      % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
alpfacd = 1; alpfacr = 1; alpfacp = 0;
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_cesqom(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,Gns,Kns,dns,ens,  tilT,tilcesq_TOT,tilom_TOT);
save([resultsdir,'/G21W_',experiment_name]); 
% negative tilcesq:
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % reverse sign of tilcesq temporarily (and also alpfacp to keep dissipation nonnegative)
nF = 2; sF = 1;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;       % amplitude of associated Legendre function of degree nF, order sF;
%tilom      = -1/2;      % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
alpfacd = 1; alpfacr = 1; alpfacp = 0;
[workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_cesqom(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,Gns,Kns,dns,ens,  tilT,tilcesq_TOT,tilom_TOT);
save([resultsdir,'/G21WN_',experiment_name]); 
tilcesq_TOT = (-1)*tilcesq_TOT; alpfacp = (-1)*alpfacp; % restore tilcesq to positive value
% %
%
end % end subfunction runcalcs_arbrot

%%
function [workAv_TOT,dissipationAv_TOT,kineticEnergyAv_TOT,potentialEnergyAv_TOT,knFsF_TOT] = solveParamSpace_cesqom(N,tilOm, alpfacd,alpfacr,alpfacp, nF,sF,Gns,Kns,dns,ens,  tilT,tilcesq_TOT,tilom_TOT)

[Nsp,Msp] = size(tilcesq_TOT);
s = sF; 
[nvec,Ntrunc] = nVec(N,s); % vector of degrees involved

parfor iitilcesq=1:Nsp, 
    percentDone = 100*iitilcesq/Nsp,
    for iitilom=1:Msp, 
        tilcesq = tilcesq_TOT(iitilcesq,iitilom);
        tilom   = tilom_TOT(iitilcesq,iitilom); 
        tilalpd = alpfacd/(tilT) * abs(tilom);
        tilalpr = alpfacr/(tilT) * abs(tilom);
        tilalpp = alpfacp/(tilT) * abs(tilom); 
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
        workAv_TOT(iitilcesq,iitilom)            = sum(calWns);
        % dissipation:
        dissipationAv_TOT(iitilcesq,iitilom)     = sum(calDns);
        % kinetic energy:
        kineticEnergyAv_TOT(iitilcesq,iitilom)   = sum(calEKns);
        % potential energy:
        potentialEnergyAv_TOT(iitilcesq,iitilom) = sum(calEPns);
        % Love number at the degree(s)/order of forcing:
        knFsF_TOT(iitilcesq,iitilom)  = knFsF;
        
    end % end iitilom loop
end % end iitilcesq loop
end % end function solveParamSpace_cesqom

end % buildParamSpace_NSR_arb




