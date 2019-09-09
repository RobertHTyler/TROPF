% This demo script shows how to calculate the tidal response (using
% tropf.m) and map variables to a global grid to make a movie of the tidal
% response over a tidal cycle (or specified time range).
%
% More specifically, it does the following: 
%
% 1) It loads a previously calculated domain of solutions and plots in Figure 1 
%    the power vs tilcesq (squared wave speed) and tilT (dissipation time
%    scale). It then brings up cross hairs onto the figure for user to
%    select a point (e.g. where a resonance peak is seen).
%
% 2) User mouse clicks on a point in the figure to obtain the parameter
%    coordinates which are assigned as tilceqChoice, tilTChoice
%    (alternatively, one can set them manually in the script by
%    uncommenting lines following "Or specify directly the coordinates:").
%
% 3) It recalculates with tropf the solution for these chosen input
%    parameters. 
%
% 4) It maps the solutions onto a global grid (2D) and further
%    calculates the fields at various times over tidal phase to create 3D
%    arrays (lon,lat,time).
%
%  5) It checks that the globe/time average of work and
%     dissipation of the gridded fields matches the same averages calculated
%     from the SH coefs. 
%
% 6) It makes movies of variables. 
% 
% 7) It converts and saves the movie to TROPF/Demos/Figs/MOV.avi
%
%
% Instructions for beginners:
%
% Step 1: First make sure you have set the root directory TROPFPATH in
% scr_startup_tropf.m; Then, set TROPFPATH in the "Startup" cell (next) to
% agree. The default here assumes TROPF is in a directory on ~/Desktop
%
% Step 2: Run this script (demo_tideMovie.m) all at once, or run by advancing
% through the cellsp.
%
% Step 3: In Figure 1, mouse cross hairs will appear. Select a point to
% obtain the fluid parameters (squared wave speed and dissipation time
% scale). Perhaps start by selecting a point in the upper right of the
% figure to see first what a near-equilibrium tide solution looks like.
% Then, for comparision, select a point within the "picket fence" of
% resonance peaks in the lower right. Open the TROPF/Figs/MOV.avi
% (e.g. with QuickTime; recommend setting "loop" in "View" tab). 
%
% Step 4: Change other input parameters and observe effects. Also change
% parameters that are plotted following instructions in the script below.



% Notes:
%
% 1) User can edit various sections of this script. This instance of the
% script is configured to load a solution domain for the retrograde
% propagating sectoral component (G22W) associated with nonsynchronous
% rapid rotation (as in Chapter Applications of manual). The solution
% domain was calculated in
% TROPF/Studies/NonSynchRotation/Limit_FastRotation and a copy moved here
% to TROPF/Demos/Data.
%
% 2) We have set N = 100 here as the default simply to speed the creation
% of the movie for this demo (tropf.m is extremely fast but the synthesizing of the movie
% frames can be slow). To gain consistency with other TROPF examples, one
% should increase this to N = 500. Note that unless a point is selected in
% the bottom right of the figure (where high-degree subharmonics are
% resonantly excited), either N = 100 or N = 500 will typically be much
% larger than needed for the accuracy required here. 
% 
% 3) If you are running Octave and this script is slow, you may try
% decreasing the value of N in the script below (e.g N = 20). Also, to
% start you could just build and plot one time frame by changing the
% specification of tiltVec (the vector of times). For example, set 
% tiltVec = 0 for one frame at time t = 0. 

% 
%
% R. Tyler, 20 April. 2019



%% % % Startup:
clear
TROPFPATH = '~/Desktop/TROPF';       % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;    % run startup script to set paths to TROPF libraries
workdir = [TROPFPATH,'/Demos'];      % working directory
cd(workdir);                         % move to working directory

fsize = 14; % font size
swCC  = 1;  % This switch is usually 1 for several function calls where the complex conjugate should be added (rather than subtracted: swCC = -1). 

 


%% Choose solution set, plot power domain, select coordinates with crosshair:

ChoiceForce = 'sectoral'  ; % example of sectoral tidal potential forcing (G22) under limit of rapid rotation
%ChoiceForce = 'tesseral' ; % example of tesseral tidal potential forcing (G22) under limit of rapid rotation

% choose solution set to load:
if ChoiceForce == 'sectoral'; load([TROPFPATH,'/Demos/Data/G22W_RD']); end
if ChoiceForce == 'tesseral'; load([TROPFPATH,'/Demos/Data/G21W_RD']); end


% plot power domain:
figvar    = 'power'; % power density
varp      = workAv_TOT; 
varpnayme = 'log$_{10}({\tilde { {\cal P }} }) $';
cl = [-10:.01:2]; cbarxticks = [-30:.5:10];
%
figure(1);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex';
cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
colormap jet; grid;
set(gca,'fontsize',fsize)
set(gca,'Xtick',10.^[-10:1:10])
caxis([-10 2]); 



% Choose a point in the power space (i.e. the parameter coordinates):
figure(1)
[tilTChoice, tilcesqChoice] = ginput(1); 
hold on; plot(tilTChoice, tilcesqChoice,'*'); hold off
%
clear *TOT varp % clear large variables or this could slow mapping below
%
%
% Or specify directly the coordinates:
% tilTChoice    = 1e1;
% tilcesqChoice = 1e-1;


%% Set other parameters and calculate solution:
% (Now that we have chosen tilcesqChoice and tilTChoice, set remaining
% input parameters and calculate tidal response with tropf.m)

% method parameters:
N     = 100; % number of SH degrees in expansion (a smaller N than the N=500 used in tropf examples can be selected to speed mapping of movie)
tilOm = 1;   % nondimensional rotation rate (should usually be 1)


% force parameters:
if ChoiceForce == 'sectoral'
nF = 2; sF = 2;s = sF; % degree and order of forcing potential 
Ntrunc = N + sF - 1;   % truncation degree
PnFsF_amp  = 3;        % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1;       % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
end
if ChoiceForce == 'tesseral'
nF = 2; sF = 1;s = sF; % degree and order of forcing potential 
Ntrunc = N + sF - 1;   % truncation degree
PnFsF_amp  = 3/2;      % amplitude of associated Legendre function of degree nF, order sF;
tilom      = -1/2;     % frequency of forcing potential  
Gns = zeros(N,1); Kns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1);
Gns((nF-sF)+1) = (1i/2) / PnFsF_amp; 
end

% fluid parameters:
tilcesq = tilcesqChoice;
tilalpd = (1/tilTChoice) * abs(tilom); 
%
tilalpr = tilalpd; tilalpp = 0;
tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;

% calculate solution using tropf.m:
[Dns,Rns,pns, calWns,calDns,calEKns,calEPnsknFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
 

%% Map solution SH coefs to grid (i.e. synthesis):

% Define grid:
DELLON = 1; DELLAT = DELLON; 
[lon,lat] = ndgrid(DELLON:DELLON:DELLON+(360-DELLON), -90:DELLAT:90); 
pixelArea = (DELLON*pi/180 * cos(lat*pi/180))*(DELLAT*pi/180);

% % Parameter settings for movie time frames:
phase0     = 0;    % phase offset (lamda_0) of obliquity
stp        = 0.05; % fractional time step
endtyme    = 1;    % fraction of orbital period we want the movie to run through
tiltVec    = [0:stp:endtyme]*2*pi; % times to sample 

% build 3D arrays of variables (lon,lat,time):
[tilD_fieldMov,  ud_fieldMov,vd_fieldMov]  = buildMovie(Dns, lon,lat,tiltVec, s,tilom,phase0,swCC);
[tilR_fieldMov,              junk1,junk2]  = buildMovie(-1i*Rns, lon,lat,tiltVec, s,tilom,phase0,swCC);
ur_fieldMov = junk2; vr_fieldMov = -junk1; clear junk1 junk2
[tilG_fieldMov]                            = buildMovie(-1i*Gns, lon,lat,tiltVec, s,tilom,phase0,swCC);
[tilp_fieldMov]                            = buildMovie(-1i*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);
[tilPTp_fieldMov]                          = buildMovie(-1i*(-1i*tilom)*pns, lon,lat,tiltVec, s,tilom,phase0,swCC);

% work field:
tilwork_fieldMov = (1/tilcesq)*tilG_fieldMov.*tilPTp_fieldMov;
% dissipation field (valid when tilalp is Rayleigh coef; otherwise need more general form):
tildissip_fieldMov = tilalpd*(ud_fieldMov.*ud_fieldMov + vd_fieldMov.*vd_fieldMov) +  tilalpr*(ur_fieldMov.*ur_fieldMov + vr_fieldMov.*vr_fieldMov) ;

%% Plot and check means of work and dissipation:

% Plot spatial distributions of time-averaged work and dissipation fields:
tilwork_field_timeAv   = mean(tilwork_fieldMov,3);
tildissip_field_timeAv = mean(tildissip_fieldMov,3);
figure(5)
contourf(lon,lat, tilwork_field_timeAv);   colorbar; colormap jet; title('time-averaged work rate density'); xlabel('longitude (degrees)'); ylabel('latitude (degrees)');
figure(6)
contourf(lon,lat, tildissip_field_timeAv); colorbar; colormap jet; title('time-averaged dissipation rate density'); xlabel('longitude (degrees)'); ylabel('latitude (degrees)');

% When these time-averaged work/dissip are further averaged over globe,
% they should match the averages obtained directly from the SH coefs. Let's
% check this:
%
% First check that the time+globe average of work and dissipation (as calculated directly from the SH coefs) match:
tilwork_timeGlobeAverage_FromSHcoefs   = sum(calWns) ; % work from SH coefs
tildissip_timeGlobeAverage_FromSHcoefs = sum(calDns) ; % dissip from SH coefs
tilwork_timeGlobeAverage_FromField   = sum(tilwork_field_timeAv(:) .* pixelArea(:)) ./ sum(pixelArea(:)) ;   % work from mapped field
tildissip_timeGlobeAverage_FromField = sum(tildissip_field_timeAv(:) .* pixelArea(:)) ./ sum(pixelArea(:)) ; % dissip from mapped field 
% the difference between average work and dissip rates obtained from SH coefs
% should be very small (and depend on number of terms N in SH expansion):
display( ['time + globe average work rate (from SH coefs) = ', num2str(tilwork_timeGlobeAverage_FromSHcoefs) ] ),
display( ['time + globe average dissipation rate (from SH coefs) = ', num2str(tildissip_timeGlobeAverage_FromSHcoefs) ] ),
display( ['maximum difference is ', num2str(abs(tilwork_timeGlobeAverage_FromSHcoefs-tildissip_timeGlobeAverage_FromSHcoefs))] ) ; 

% Now check that the time+globe averages of the gridded work/dissip fields match the
% time+globe averages calculated directly from the SH coefs:
%
% The differences should be small (but the accuracy of the mapped fields
% depends on the grid resolution, so the differences may not be as small as
% the differences using only SH coefs---which are generally more accurate):
display( ['time + globe average work rate (from field) = ', num2str(tilwork_timeGlobeAverage_FromField) ] ),
display( ['time + globe average dissipation rate (from field) = ', num2str(tildissip_timeGlobeAverage_FromField) ] ),
display( ['maximum difference in work rate (field vs SH) = ', num2str(abs(tilwork_timeGlobeAverage_FromField-tilwork_timeGlobeAverage_FromSHcoefs))] ) ; 
display( ['maximum difference in dissip rate (field vs SH) = ', num2str(abs(tildissip_timeGlobeAverage_FromField-tildissip_timeGlobeAverage_FromSHcoefs))] ) ; 


 


%% Make a movie:

% Choose variables for movie:
%
% choose variable to contour:
%varpMov   = tilD_fieldMov;  % divergent flow potential
%varpMov   = tilR_fieldMov;  % rotational flow potential
varpMov   = tilp_fieldMov;   % dynamic pressure variable (think surface displacement in simplest barotropic case)
%varpMov   = tilwork_fieldMov;   % work rate density
%varpMov   = tildissip_fieldMov; % dissipation rate density 
%
% choose variables to quiver:
varpMov_E = ud_fieldMov + ur_fieldMov; % "eastward" (prograde) velocity
varpMov_N = vd_fieldMov + vr_fieldMov; % "northward" velocity
%varpMov_E = (ud_fieldMov + ur_fieldMov).*tilp_fieldMov; % eastward wave momentum
%varpMov_N = (vd_fieldMov + vr_fieldMov).*tilp_fieldMov; % northward wave momentum


% choose name of contoured variable for title:
varpnayme = '$ {\tilde { p } } $';



% make movie: 
[Mm,Nn,Tt] = size(varpMov);
for ii = 1:Tt
varp   = squeeze(varpMov(:,:,ii));
varp_E = squeeze(varpMov_E(:,:,ii));
varp_N = squeeze(varpMov_N(:,:,ii));
fig10 = figure(10);fig10.Units = 'inches';
fig10.OuterPosition = [8,5,12,8];clf
colormap(othercolor('BuDRd_18'));
cax = [-1 1]*max(abs(varp(:)));
contour3(lon,lat,varp,[-1:0.02:1]*cax(1)); xlabel('longitude (deg.)'); ylabel('latitude (deg.)');
caxis(cax); zlims = cax*3; xlims = [lon(1) lon(end)]; ylims = [lat(1) lat(end)];
set(gca,'xlim',xlims,'ylim',ylims,'zlim',zlims,'Box','off'); view(-15,40);
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
hold on;
stplon = 10; stplat = stplon; 
lonjunk = lon(1:stplon:end, 1:stplat:end);
latjunk = lat(1:stplon:end, 1:stplat:end);
uvarpjunk = varp_E(1:stplon:end, 1:stplat:end);
vvarpjunk = varp_N(1:stplon:end, 1:stplat:end);
quiver3(lonjunk,latjunk,0*latjunk+zlims(1), uvarpjunk,vvarpjunk,0*latjunk)
title([varpnayme, ' and flow vectors'],'Interpreter','latex')
hold off
MOV(ii) = getframe(fig10);
end


% save movie
vidjunk = VideoWriter('./Figs/MOV.avi','Uncompressed AVI'); vidjunk.FrameRate = 5;
open(vidjunk); writeVideo(vidjunk,MOV); close(vidjunk);



