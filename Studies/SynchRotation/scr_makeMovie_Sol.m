% Example script does the following:
% 1) Loads a previously calculated solution (power) space and makes contour plot
% 2) User mouse clicks on a point in the figure to obtain the parameter coordinates
% 3) Recalculates with tropf the solution for these chosen input parameters
% 4) Maps the solutions onto a grid (2D) and further calculates at various times over tidal phase to create 3D arrays (lon,lat,time)
% 5) Checks that the globe/time average of work and dissipation of the gridded fields matches the same averages calculated from the SH coefs
% 6) Makes moves of variables
% 7) Converts and saves an avi movie
% Notes:
% User can edit various sections. This instance of this script is configured to load solutions from the NonSynchRotation directory

% 
%
% R. Tyler, 10 Jan. 2019



%% % % build movie of tidal potential  
clear
TROPFPATH = '~/Desktop/TROPF';       % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;    % run startup script to set paths to TROPF libraries
workdir = [TROPFPATH,'/Studies/SynchRotation'];   % working directory
cd(workdir);                         % move to working directory



fsize = 14; % font size for plots
swCC  = 1;  % This is usually 1 for several function calls where the complex conjugate should be added (rather than subtracted: swCC = -1). 

 
%% Choose a power space figure and ginput to get coordinates of a point:


% Choose a point in the power space (i.e. the parameter coordinates):
if 1==2
figure(101)
[tilTChoice, tilcesqChoice] = ginput(1); 
hold on; plot(tilTChoice, tilcesqChoice,'*'); hold off
end

% or specify directly the coordinates:
if 1 == 1
tilTChoice    = 1e2;
tilcesqChoice = 1e-2;
end


clear *TOT varp % clear large variables or this will slow mapping below



%% Calculate solutions (SH coefs) and map fields onto lon,lat,time:
%
% Set method parameters common to all components:
N     = 10;
tilOm = 1; 
%
% Define grid:
DELLON = 1; DELLAT = DELLON; 
[lon,lat] = ndgrid(DELLON:DELLON:DELLON+(360-DELLON), -90:DELLAT:90); 
pixelArea = (DELLON*pi/180 * cos(lat*pi/180))*(DELLAT*pi/180);
%
% % Parameter settings for movie frames:
phase0     = 0;    % phase offset (lamda_0) of obliquity
stp        = 0.05; % fractional time step
endtyme    = 1;    % fraction of orbital period we want the movie to run through
tiltVec    = [0:stp:endtyme]*2*pi *(2); % times to sample ( *(2) so that the tilom=1/2 components run through full cycle)
% 
% %% Calculate solutions (SH coefs) and map fields onto lon,lat,time:
% %
% % Set method parameters common to all components:
% N     = 20;  % number of terms in SH expansion
% tilOm = 1;    % nondimensional rotation rate
% %
% % Set fluid parameters common to all components:
% tilcesq      = tilcesqChoice; % squared nondimensional wavespeed
% abstilom     = 1/2;           % nondimensional amplitude of tilom for all synch rotation components is 1/2 (sign, however, differs)
% tilalpd      = (1/tilTChoice)*abstilom;   % nondimensional dissipation coef for divergent flow
% tilalpr      = tilalpd;                   % nondimensional dissipation coef for rotational flow
% tilalpp      = 0;                         % nondimensional dissipation coef for 
% tilnusqnsPos =  (1 + 1i*tilalpp/abstilom)./tilcesq;      % (for positive tilom); nondimensional complex slowness parameter vector; note that tilcesq and tilalpp enter tropf only through the slowness parameter.
% tilnusqnsNeg =  (1 + 1i*tilalpp/(-1*abstilom))./tilcesq; % (for negative tilom); nondimensional complex slowness parameter vector; note that tilcesq and tilalpp enter tropf only through the slowness parameter. This is also the cc of tilnusqnsPos.
% 
% % Define grid:
% DELLON = 1; DELLAT = DELLON; 
% [lon,lat] = ndgrid(DELLON:DELLON:DELLON+(360-DELLON), -90:DELLAT:90); 
% pixelArea = (DELLON*pi/180 * cos(lat*pi/180))*(DELLAT*pi/180);
% %
% % Set parameters for time expansion and movie frames:
% phase0     = 0;    % phase offset (lamda_0) of obliquity
% stp        = 0.05; % fractional time step
% endtyme    = 1;    % fraction of period we want the movie to run through
% tiltVec    = [0:stp:endtyme]*2*pi *(1/abstilom); % vector of times at which to sample
% 
 
% Eccentricity libration G20:
%
% tropf solve for SH coefs:
nF = 2; sF = 0;s=sF;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); tilalpr = tilalpd; tilalpp = 0;        
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;
% Rename some solution variables:
PnFsF_amp20=PnFsF_amp;Dns20=Dns;Rns20=Rns;pns20=pns; calWns20=calWns;calDns20=calDns;calEKns20=calEKns;calEPns20=calEPns;knFsF20=knFsF;
tilD_fieldMov20=tilD_fieldMov; ud_fieldMov20=ud_fieldMov;vd_fieldMov20=vd_fieldMov;
tilR_fieldMov20=tilR_fieldMov; ur_fieldMov20=ur_fieldMov;vr_fieldMov20=vr_fieldMov;
tilG_fieldMov20=tilG_fieldMov; tilp_fieldMov20=tilp_fieldMov;
tilwork_fieldMov20=tilwork_fieldMov;
tilp_fieldMov20=tilp_fieldMov;tildissip_fieldMov20=tildissip_fieldMov;
%
%
% Eccentricity libration (westward prop component) G22W:
%
% tropf solve for SH coefs:
nF = 2; sF = 2;s=sF       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); tilalpr = tilalpd; tilalpp = 0;
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;
% Rename some solution variables:
PnFsF_amp22W=PnFsF_amp;Dns22W=Dns;Rns22W=Rns;pns22W=pns; calWns22W=calWns;calDns22W=calDns;calEKns22W=calEKns;calEPns22W=calEPns;knFsF22W=knFsF;
tilD_fieldMov22W=tilD_fieldMov; ud_fieldMov22W=ud_fieldMov;vd_fieldMov22W=vd_fieldMov;
tilR_fieldMov22W=tilR_fieldMov; ur_fieldMov22W=ur_fieldMov;vr_fieldMov22W=vr_fieldMov;
tilG_fieldMov22W=tilG_fieldMov; tilp_fieldMov22W=tilp_fieldMov;
tilwork_fieldMov22W=tilwork_fieldMov;
tilp_fieldMov22W=tilp_fieldMov;tildissip_fieldMov22W=tildissip_fieldMov;
%
%
% Eccentricity libration (westward prop component) G22E:
%
% tropf solve for SH coefs:
nF = 2; sF = 2;s=sF       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3;       % amplitude of associated Legendre function of degree nF, order sF;
tilom      = 1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); tilalpr = tilalpd; tilalpp = 0;
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;
% Rename some solution variables:
PnFsF_amp22E=PnFsF_amp;Dns22E=Dns;Rns22E=Rns;pns22E=pns; calWns22E=calWns;calDns22E=calDns;calEKns22E=calEKns;calEPns22E=calEPns;knFsF22E=knFsF;
tilD_fieldMov22E=tilD_fieldMov; ud_fieldMov22E=ud_fieldMov;vd_fieldMov22E=vd_fieldMov;
tilR_fieldMov22E=tilR_fieldMov; ur_fieldMov22E=ur_fieldMov;vr_fieldMov22E=vr_fieldMov;
tilG_fieldMov22E=tilG_fieldMov; tilp_fieldMov22E=tilp_fieldMov;
tilwork_fieldMov22E=tilwork_fieldMov;
tilp_fieldMov22E=tilp_fieldMov;tildissip_fieldMov22E=tildissip_fieldMov;
%
%
% Obliquity libration (westward prop component) G21W:
%
% positive tilcesq:
nF = 2; sF = 1;s=sF;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); tilalpr = tilalpd; tilalpp = 0;
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;
% Rename some solution variables:
PnFsF_amp21W=PnFsF_amp;Dns21W=Dns;Rns21W=Rns;pns21W=pns; calWns21W=calWns;calDns21W=calDns;calEKns21W=calEKns;calEPns21W=calEPns;knFsF21W=knFsF;
tilD_fieldMov21W=tilD_fieldMov; ud_fieldMov21W=ud_fieldMov;vd_fieldMov21W=vd_fieldMov;
tilR_fieldMov21W=tilR_fieldMov; ur_fieldMov21W=ur_fieldMov;vr_fieldMov21W=vr_fieldMov;
tilG_fieldMov21W=tilG_fieldMov; tilp_fieldMov21W=tilp_fieldMov;
tilwork_fieldMov21W=tilwork_fieldMov;
tilp_fieldMov21W=tilp_fieldMov;tildissip_fieldMov21W=tildissip_fieldMov;
%
%
% Obliquity libration (westward prop component) G21E:
%
% positive tilcesq:
nF = 2; sF = 1;s=sF;       % degree and order of forcing potential 
Ntrunc = N + sF - 1;  % truncation degree
PnFsF_amp  = 3/2;     % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;    % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); tilalpr = tilalpd; tilalpp = 0;
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);;
% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;
% Rename some solution variables:
PnFsF_amp21E=PnFsF_amp;Dns21E=Dns;Rns21E=Rns;pns21E=pns; calWns21E=calWns;calDns21E=calDns;calEKns21E=calEKns;calEPns21E=calEPns;knFsF21E=knFsF;
tilD_fieldMov21E=tilD_fieldMov; ud_fieldMov21E=ud_fieldMov;vd_fieldMov21E=vd_fieldMov;
tilR_fieldMov21E=tilR_fieldMov; ur_fieldMov21E=ur_fieldMov;vr_fieldMov21E=vr_fieldMov;
tilG_fieldMov21E=tilG_fieldMov; tilp_fieldMov21E=tilp_fieldMov;
tilwork_fieldMov21E=tilwork_fieldMov;
tilp_fieldMov21E=tilp_fieldMov;tildissip_fieldMov21E=tildissip_fieldMov;
 






%% Some checks to see that integrals of mapped fields and SH coefs agree:
% Plot distributions of work and dissipation; Also note that full means should be close (there's a difference due to pixel weighting dependence with lat)
figure(44)
junk1 = pixelArea.*mean(tilwork_fieldMov,3); 
junk2 = pixelArea.*mean(tildissip_fieldMov,3); 
% check that globe/time average work and dissipation from either SH coefs or mapped fields match:
[sum(calWns),sum(junk1(:))./sum(pixelArea(:))]
[sum(calDns),sum(junk2(:))./sum(pixelArea(:))]

figure(44)
contourf(lon,lat, junk1); colorbar
figure(45)
contourf(lon,lat, junk2); colorbar

return
%% Compile total components, scaling for relative amplitudes:
ecc = 9e-3;   % eccentricity parameter
obl = 1.7e-3;   % obliquity parameter
% factors in components of synch grav potential (times -1):
fac21W = 1/2*obl; fac21E = 1/2*obl;                     
fac20 = -3/2*ecc; fac22W = -1/8*ecc; fac22E = 7/8*ecc;  

% Override above values to turn off some components;
%fac20=0;
%fac22W=1;
%fac22E=1;


tilG_fieldMov = fac20*PnFsF_amp20*tilG_fieldMov20 + fac22W*PnFsF_amp22W*tilG_fieldMov22W + fac22E*PnFsF_amp22E*tilG_fieldMov22E + fac21W*PnFsF_amp21W*tilG_fieldMov21W + fac21E*PnFsF_amp21E*tilG_fieldMov21E;
tilD_fieldMov = fac20*PnFsF_amp20*tilD_fieldMov20 + fac22W*PnFsF_amp22W*tilD_fieldMov22W + fac22E*PnFsF_amp22E*tilD_fieldMov22E + fac21W*PnFsF_amp21W*tilD_fieldMov21W + fac21E*PnFsF_amp21E*tilD_fieldMov21E;
tilR_fieldMov = fac20*PnFsF_amp20*tilR_fieldMov20 + fac22W*PnFsF_amp22W*tilR_fieldMov22W + fac22E*PnFsF_amp22E*tilR_fieldMov22E + fac21W*PnFsF_amp21W*tilR_fieldMov21W + fac21E*PnFsF_amp21E*tilR_fieldMov21E;
tilp_fieldMov = fac20*PnFsF_amp20*tilp_fieldMov20 + fac22W*PnFsF_amp22W*tilp_fieldMov22W + fac22E*PnFsF_amp22E*tilp_fieldMov22E + fac21W*PnFsF_amp21W*tilp_fieldMov21W + fac21E*PnFsF_amp21E*tilp_fieldMov21E;
tilwork_fieldMov = fac20*PnFsF_amp20*tilwork_fieldMov20 + fac22W*PnFsF_amp22W*tilwork_fieldMov22W + fac22E*PnFsF_amp22E*tilwork_fieldMov22E + fac21W*PnFsF_amp21W*tilwork_fieldMov21W + fac21E*PnFsF_amp21E*tilwork_fieldMov21E;
tildissip_fieldMov = fac20*PnFsF_amp20*tildissip_fieldMov20 + fac22W*PnFsF_amp22W*tildissip_fieldMov22W + fac22E*PnFsF_amp22E*tildissip_fieldMov22E + fac21W*PnFsF_amp21W*tildissip_fieldMov21W + fac21E*PnFsF_amp21E*tildissip_fieldMov21E;
ud_fieldMov = fac20*PnFsF_amp20*ud_fieldMov20 + fac22W*PnFsF_amp22W*ud_fieldMov22W + fac22E*PnFsF_amp22E*ud_fieldMov22E + fac21W*PnFsF_amp21W*ud_fieldMov21W + fac21E*PnFsF_amp21E*ud_fieldMov21E;
vd_fieldMov = fac20*PnFsF_amp20*vd_fieldMov20 + fac22W*PnFsF_amp22W*vd_fieldMov22W + fac22E*PnFsF_amp22E*vd_fieldMov22E + fac21W*PnFsF_amp21W*vd_fieldMov21W + fac21E*PnFsF_amp21E*vd_fieldMov21E;
ur_fieldMov = fac20*PnFsF_amp20*ur_fieldMov20 + fac22W*PnFsF_amp22W*ur_fieldMov22W + fac22E*PnFsF_amp22E*ur_fieldMov22E + fac21W*PnFsF_amp21W*ur_fieldMov21W + fac21E*PnFsF_amp21E*ur_fieldMov21E;
vr_fieldMov = fac20*PnFsF_amp20*vr_fieldMov20 + fac22W*PnFsF_amp22W*vr_fieldMov22W + fac22E*PnFsF_amp22E*vr_fieldMov22E + fac21W*PnFsF_amp21W*vr_fieldMov21W + fac21E*PnFsF_amp21E*vr_fieldMov21E;
 




 
%% Make a movie:
% choose variable to contour3:
%varpMov   = tilG_fieldMov;
%varpMov   = tilR_fieldMov;
varpMov   = tilp_fieldMov;
%varpMov   = tilwork_fieldMov;
%varpMov   = tildissip_fieldMov;
%varpMov   = tilp_fieldMov.*varpMov_E;
% choose variables to quiver:
varpMov_E = ud_fieldMov + ur_fieldMov; % "eastward" (prograde) velocity
varpMov_N = vd_fieldMov + vr_fieldMov; % "northward" velocity
%varpMov_E = (ud_fieldMov + ur_fieldMov).*tilp_fieldMov; % eastward wave momentum
%varpMov_N = (vd_fieldMov + vr_fieldMov).*tilp_fieldMov; % northward wave momentum



[Mm,Nn,Tt] = size(varpMov);

for ii = 1:Tt
varp   = squeeze(varpMov(:,:,ii));
varp_E = squeeze(varpMov_E(:,:,ii));
varp_N = squeeze(varpMov_N(:,:,ii));
 
%fig10 = figure(10);
fig10 = figure(10);fig10.Units = 'inches';
%fig10.OuterPosition = [12,2,8,10];clf
fig10.OuterPosition = [8,5,12,8];clf
%fig10.OuterPosition = [12,2,8,10];clf
colormap(othercolor('BuDRd_18'));

cax = [-1 1]*max(abs(varpMov(:)))
contour3(lon,lat,varp,[-1:0.02:1]*cax(1)); xlabel('longitude (deg.)'); ylabel('latitude (deg.)');
caxis(cax); zlims = cax*3; xlims = [lon(1) lon(end)]; ylims = [lat(1) lat(end)];
set(gca,'xlim',xlims,'ylim',ylims,'zlim',zlims,'Box','off'); view(-15,25);
title('varpMov'); junk = colorbar; junk.Label.String = 'varpMov units';

hold on;
stplon = 10; stplat = stplon; 
lonjunk = lon(1:stplon:end, 1:stplat:end);
latjunk = lat(1:stplon:end, 1:stplat:end);
uvarpjunk = varp_E(1:stplon:end, 1:stplat:end);
vvarpjunk = varp_N(1:stplon:end, 1:stplat:end);
quiver3(lonjunk,latjunk,0*latjunk+zlims(1), uvarpjunk,vvarpjunk,0*latjunk)
%scjunk = 20/max(abs(uvarpjunk(:)));
%quiver3(lonjunk,latjunk,0*latjunk+zlims(1), scjunk*uvarpjunk,scjunk*vvarpjunk,0*latjunk,'Autoscale','off')
hold off

MOV(ii) = getframe(fig10)
end


%%

% save movie
savemovie = 1
if savemovie
vidjunk = VideoWriter('./Figs/MOV.avi','Uncompressed AVI'); vidjunk.FrameRate = 5;
open(vidjunk); writeVideo(vidjunk,MOV); close(vidjunk);
end

%% Plot time means of some variables:


%%
%figvar    = 'varp'; % power density
%varp      = mean(tilwork_fieldMov,3); 
varp      = mean(tildissip_fieldMov,3); 
%varp      = mean(abs(tilp_fieldMov),3); 
%varpnayme = 'log$_{10}({\tilde { {\cal P }} }) $';
%cl = [-10:.1:15]; cbarxticks = [-2:1:15];
%cl = [-10:.01:2]; cbarxticks = [-30:.5:10];
%
figure(100);clf
contourf(lon,lat,varp,'linestyle','none'); shading flat; 
xlabel('longitude (deg.)'); ylabel('latitude (deg.)'); 
cbar = colorbar('Eastoutside');%,'xtick',cbarxticks); cbar.Label.Interpreter = 'latex';
%cbar.Label.String = varpnayme;
set(gca,'fontsize',fsize)
%set(gca,'Xtick',10.^[-10:1:10])
%caxis([-10 2]); 
%scr_reduceWhiteMargins
%print([workdir,'/Figs/',figvar,'_',component],'-dpng');


