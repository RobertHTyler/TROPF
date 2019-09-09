% This demo script uses tropf.m to calculate the tidal power and response (e.g.
% 'k_2'), and determine if it is close to the equilibrium-tide solution.


% Instructions for beginners:
%
% Step 1: First make sure you have set the root directory TROPFPATH in
% scr_startup_tropf.m; Then, set TROPFPATH in the "Startup" cell (next) to
% agree. The default here assumes TROPF is in a directory on ~/Desktop
%
% Step 2: Run this script (demo_power.m) all at once, or run by advancing
% through the cells to follow the execution better.
%
% Step 3: See notes in the "Inspect Solution" cells below. You should first
% see that the average tidal work rate density is equal to the average
% tidal dissipation rate. Either can then be regarded as the "tidal power
% density". To dimensionalize solutions (e.g. to obtain power in W/m^3) for
% a specific application, see dimensionalization factors in Table 1 of the
% TROPF manual. Then look at the admittance (will be near 1 for an
% equilibrium-tide response). Then look at the figures with the
% distributions with respect to degree.
%
% Step 4: Edit input parameters to see how tidal response changes.
% Suggestions:
%
% a) Change force parameters: This can be done by uncommenting a different
% "ChoiceForce" for the two examples provided (which are the sectoral and
% tesseral components of the force on a rapid-rotation nonsychronously
% rotating fluid, as in the Applications Chapter (Non-Synchronous Rapid
% Rotation section) of the TROPF manual). Alternatively, one may hand edit
% the "Force parameters" paragraph to reflect an arbitrary choice for the
% force.
%
% b) Change fluid parameters: For example set the squared wave speed
% (tilcesq) to tilcesq < 1 for a non-equilibrium response. For a
% non-equilbrium response, you may also or instead increase the (here,
% Rayleigh-form) dissipation by setting tilalpd = tilalpr > 1. To see
% resonantly forced scenario, assign tilceq to be one of the eigenvalues
% described in the Applications Chapter (Non-Synchronous Rapid Rotation
% section) of the manual (for sectoral forcing, for example, you could
% assign the first eigenvalue of the symmetric (even order) mode as
% follows: tilcesq = 9.98e-2


% Notes:
%
% 1) All parameters, variables, and operators have been non-dimensionalized
% following TROPF manual.
%
% 2) Abbreviations: vec (vec); SH (spherical harmonic); coef (coefficient).
%
% 3) In the equilibrium-tide solution, the wave speeds in the fluid are
% fast enough to quasi-statically adjust the fluid to the tidal forcing
% potential. This also requires that dissipation be weak. 
%
% 4) See further notes in cells below.
% %

% R. Tyler, 10 Feb, 2019


%% % % Startup:

clear
TROPFPATH  = '~/Desktop/TROPF';       % set path of TROPF root directory (should agree with TROPFPATH set in scr_startup_tropf.m)
cd(TROPFPATH); scr_startup_tropf;     % run startup script to set paths to TROPF libraries
workdir    = [TROPFPATH,'/Demos'];    % working directory
cd(workdir)



%% % % Set input parameters and solve using tropf.m:

% Method parameters:
tilOm = 1   ; % nondim rotation rate (\Omega / \Omega_s); usually set = 1
N     = 500 ; % number of terms in SH expansion



% Force parameters:
%
% Choose the sectoral or tesseral tidal force component:
ChoiceForce = 'sectoral';
%ChoiceForce = 'tesseral';
%
% Sectoral tidal force potential (G22W):
if ChoiceForce == 'sectoral'
nF = 2; sF = 2;       % degree and order of forcing potential 
s  = sF;              % the order (s) in the expansion is the same order (sF) as the forcing
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1;      % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1); % initialize all force coefs to zero
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; % sets unit amplitude sectoral tidal force potential 
end
% Tesseral tidal force potential (G21W):
if ChoiceForce == 'tesseral'
nF = 2; sF = 1;       % degree and order of forcing potential 
s  = sF;              % the order (s) in the expansion is the same order (sF) as the forcing
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1); % initialize all force coefs to zero
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; % sets unit amplitude sectoral tidal force potential 
end



% Fluid parameters:
tilalpd = 0.01; tilalpr = 0.01; tilalpp = 0;  % dissipation coefs 
tilcesq = 1e4;                                % squared wave speed
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq; % squared slowness



% Solve for the SH coefs of the tidal response:
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,sF,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);


%% % % Inspect solution:

% tidal power:
averageWork = sum(calWns)   % average work rate density
averageDiss = sum(calDns)   % average dissipation rate density
maxDiff     = max(abs(averageWork(:) - averageDiss(:)))  % the difference should be very small
tidalPower  = averageWork   % provided maxDiff is small, take the tidal power to be equal to either the average work or dissipation rate density


% tidal response admittance at degree and order of forcing (sometimes regarded as complex 'k_2'):
knFsF  % will be close to 1 for an equlibrium tide response; will have large imaginary part if strong dissipation


%% % % Inspect solution (cont.):
%
% While an equilibrium tide requires knFsF = 1, this alone is not a sufficient
% condition to claim an equilibrium-tide response as there may also be
% tidal response energy at degrees other than that of the forcing. To look
% at this, let's plot below the distributions with respect to degree. If
% the pressure variable shows amplitude only at the same degree as the
% forcing, then we can claim the equilibrium-tide response. 


% % Get vec of SH degrees (n = s, s+1, s+2 ... Ntrunc):
[nvec,Ntrunc] = nVec(N,s);

figure(1); subplot(311)
bar(nvec(1:10), abs(Dns(1:10))); xlabel('degree'); ylabel('divergent flow potential')
figure(1); subplot(312)
bar(nvec(1:10), abs(Rns(1:10))); xlabel('degree'); ylabel('rotational flow potential')
figure(1); subplot(313)
bar(nvec(1:10), abs(pns(1:10))); xlabel('degree'); ylabel('dynamic pressure')

figure(2); subplot(211)
bar(nvec(1:10), calWns(1:10)); xlabel('degree'); ylabel('work')
figure(2); subplot(212)
bar(nvec(1:10), calDns(1:10)); xlabel('degree'); ylabel('dissipation')

