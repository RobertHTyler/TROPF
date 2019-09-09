% Plot power domain (i.e. plot power as a function of tilT, tilcesq):
% Assumes the solutions have first been built (e.g. using
% buildParamSpace_*.m)
%
%
% Notes: 
%
% 1) Here the variable "component" refers to the individual complete
% solution data set we wish to plot (which may vary with dissipation scheme
% assumptions for example). But "component" is also used in the comments to
% refer to the component of forcing (because the synchronous rotation case
% needs to consider eccentricity and/or obliquity driven tides which
% involve combining multiple solutions (e.g. eccentricity involves three
% components (G20, G22W, G22E, which are respectively the radial; n=2,s=2 westward propagating;
% and n=2,s=2 eastward propagating components). An N is added (e.g. G20N)
% to indicate that the data set is for negative tilcesq. 

%% Startup:

clear
TROPFPATH  = '~/Desktop/TROPF';     % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;   % run startup script to set paths to TROPF libraries
resultsdir = [TROPFPATH,'/Studies/SynchRotation/Results'];  % working directory
workdir    = [TROPFPATH,'/Studies/SynchRotation'];  % working directory
cd(workdir);                        % move to working directory

fsize = 14; % font size for plots
 
 
%% Choose solution set to load

%component = 'RD'; % Rayleigh dissipation assumption
component = 'NR'; % Newtonian radiation assumption

% switch negtilcesq_sw set to 1 if we should load the results for negative
% tilcesq; otherwise set to 0 for the positive tilcesq results:
for negtilcesq_sw = [0 1] % loop to plot both positive and negative tilcesq data sets
%close all

% Consider the scaled combination of all three components (G20, G22W, G22E)
% of the eccentricty tides:
ecc = 1;   % eccentricity parameter (needed for scaling)
obl = 1;   % obliquity parameter (radians; needed for scaling)
% factors in components of synch grav potential (times -1):
fac21W = 1/2*obl; fac21E = 1/2*obl;                     
fac20 = -3/2*ecc; fac22W = -1/8*ecc; fac22E = 7/8*ecc;  



if negtilcesq_sw == 0
% Load eccentricity components:
load([resultsdir,'/','G20_',component]);
workG20            = (fac20*PnFsF_amp)^2*workAv_TOT;
dissipationG20     = (fac20*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG20   = (fac20*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG20 = (fac20*PnFsF_amp)^2*potentialEnergyAv_TOT;
 
load([resultsdir,'/','G22W_',component]);
workG22W            = (fac22W*PnFsF_amp)^2*workAv_TOT;
dissipationG22W     = (fac22W*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG22W   = (fac22W*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG22W = (fac22W*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G22E_',component]);
workG22E            = (fac22E*PnFsF_amp)^2*workAv_TOT;
dissipationG22E     = (fac22E*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG22E   = (fac22E*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG22E = (fac22E*PnFsF_amp)^2*potentialEnergyAv_TOT;
 
% Load obliquity components:
load([resultsdir,'/','G20_',component]);
workG20            = (fac20*PnFsF_amp)^2*workAv_TOT;
dissipationG20     = (fac20*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG20   = (fac20*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG20 = (fac20*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G21W_',component]);
workG21W            = (fac21W*PnFsF_amp)^2*workAv_TOT;
dissipationG21W     = (fac21W*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG21W   = (fac21W*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG21W = (fac21W*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G21E_',component]);
workG21E            = (fac21E*PnFsF_amp)^2*workAv_TOT;
dissipationG21E     = (fac21E*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG21E   = (fac21E*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG21E = (fac21E*PnFsF_amp)^2*potentialEnergyAv_TOT;

elseif negtilcesq_sw == 1
% Load eccentricity components:
load([resultsdir,'/','G20N_',component]);
workG20            = (fac20*PnFsF_amp)^2*workAv_TOT;
dissipationG20     = (fac20*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG20   = (fac20*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG20 = (fac20*PnFsF_amp)^2*potentialEnergyAv_TOT;
 
load([resultsdir,'/','G22WN_',component]);
workG22W            = (fac22W*PnFsF_amp)^2*workAv_TOT;
dissipationG22W     = (fac22W*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG22W   = (fac22W*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG22W = (fac22W*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G22EN_',component]);
workG22E            = (fac22E*PnFsF_amp)^2*workAv_TOT;
dissipationG22E     = (fac22E*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG22E   = (fac22E*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG22E = (fac22E*PnFsF_amp)^2*potentialEnergyAv_TOT;
 
% Load obliquity components:
load([resultsdir,'/','G20N_',component]);
workG20            = (fac20*PnFsF_amp)^2*workAv_TOT;
dissipationG20     = (fac20*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG20   = (fac20*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG20 = (fac20*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G21WN_',component]);
workG21W            = (fac21W*PnFsF_amp)^2*workAv_TOT;
dissipationG21W     = (fac21W*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG21W   = (fac21W*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG21W = (fac21W*PnFsF_amp)^2*potentialEnergyAv_TOT;

load([resultsdir,'/','G21EN_',component]);
workG21E            = (fac21E*PnFsF_amp)^2*workAv_TOT;
dissipationG21E     = (fac21E*PnFsF_amp)^2*dissipationAv_TOT;
kineticEnergyG21E   = (fac21E*PnFsF_amp)^2*kineticEnergyAv_TOT;
potentialEnergyG21E = (fac21E*PnFsF_amp)^2*potentialEnergyAv_TOT;
end








% Choose total to be just the  G20 eccentricity components:
if 1 == 2
workAv_TOT            = workG20;
dissipationAv_TOT     = dissipationG20;
kineticEnergyAv_TOT   = kineticEnergyG20;
potentialEnergyAv_TOT = potentialEnergyG20;
if negtilcesq_sw == 0, compsub = 'G20', elseif negtilcesq_sw == 1 compsub = 'G20N', end
end



% Choose total to be just the G22W eccentricity components:
if 1 == 2
workAv_TOT            = workG22W;
dissipationAv_TOT     = dissipationG22W;
kineticEnergyAv_TOT   = kineticEnergyG22W;
potentialEnergyAv_TOT = potentialEnergyG22W;
if negtilcesq_sw == 0, compsub = 'G22W', elseif negtilcesq_sw == 1 compsub = 'G22WN', end
end



% Choose total to be just the G22E eccentricity components:
if 1 == 2
workAv_TOT            = workG22E;
dissipationAv_TOT     = dissipationG22E;
kineticEnergyAv_TOT   = kineticEnergyG22E;
potentialEnergyAv_TOT = potentialEnergyG22E;
if negtilcesq_sw == 0, compsub = 'G22E', elseif negtilcesq_sw == 1 compsub = 'G22EN', end
end



% Choose total to be just the G21W obliquity components:
if 1 == 2
workAv_TOT            = workG21W;
dissipationAv_TOT     = dissipationG21W;
kineticEnergyAv_TOT   = kineticEnergyG21W;
potentialEnergyAv_TOT = potentialEnergyG21W;
if negtilcesq_sw == 0, compsub = 'G21W', elseif negtilcesq_sw == 1 compsub = 'G21WN', end
end
 


% Choose total to be just the G21E obliquity components:
if 1 == 2
workAv_TOT            = workG21E;
dissipationAv_TOT     = dissipationG21E;
kineticEnergyAv_TOT   = kineticEnergyG21E;
potentialEnergyAv_TOT = potentialEnergyG21E;
if negtilcesq_sw == 0, compsub = 'G21E', elseif negtilcesq_sw == 1 compsub = 'G21EN', end
end
 

 




% Choose total to be just the eccentricity components:
if 1 == 2
workAv_TOT            = workG20 + workG22W + workG22E;
dissipationAv_TOT     = dissipationG20 + dissipationG22W + dissipationG22E;
kineticEnergyAv_TOT   = kineticEnergyG20 + kineticEnergyG22W + kineticEnergyG22E;
potentialEnergyAv_TOT = potentialEnergyG20 + potentialEnergyG22W + potentialEnergyG22E;
if negtilcesq_sw == 0, compsub = 'Ecc', elseif negtilcesq_sw == 1 compsub = 'EccN', end
end


% Choose total to be just the obliquity components:
if 1 == 1
workAv_TOT            = workG21W + workG21E;
dissipationAv_TOT     = dissipationG21W + dissipationG21E;
kineticEnergyAv_TOT   = kineticEnergyG21W + kineticEnergyG21E;
potentialEnergyAv_TOT = potentialEnergyG21W + potentialEnergyG21E;
if negtilcesq_sw == 0, compsub = 'Obl', elseif negtilcesq_sw == 1 compsub = 'OblN', end
end
 


% Choose total to be just the eccentricity+obliquity components:
if 1 == 2
workAv_TOT            = workG20 + workG22W + workG22E  +  workG21W + workG21E;
dissipationAv_TOT     = dissipationG20 + dissipationG22W + dissipationG22E  +  dissipationG21W + dissipationG21E;
kineticEnergyAv_TOT   = kineticEnergyG20 + kineticEnergyG22W + kineticEnergyG22E  +  kineticEnergyG21W + kineticEnergyG21E;
potentialEnergyAv_TOT = potentialEnergyG20 + potentialEnergyG22W + potentialEnergyG22E  +  potentialEnergyG21W + potentialEnergyG21E;
if negtilcesq_sw == 0, compsub = 'EccObl', elseif negtilcesq_sw == 1 compsub = 'EccOblN', end
end


%% Check some equivalences:

mabs(workAv_TOT,dissipationAv_TOT) % should be small since average work and dissipation are equal

%mabs(dissipationAv_TOT,  2*tilalpd*kineticEnergyAv_TOT)



%% Plots:

 

%%
figvar    = 'power'; % power density
varp      = workAv_TOT; 
varpnayme = 'log$_{10}({\tilde { {\cal P }} }) $';
cl = [-10:.01:2]; cbarxticks = [-30:.5:10];
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');


%%
figvar    = 'tilcesqpower'; % power density multiplied by tilcesq (prop. to layer thickness)
varp      = abs(tilcesq_TOT).*workAv_TOT; 
varpnayme = 'log$_{10}(|{\tilde c_e}^2|{\tilde { {\cal P }} }) $';
cl = [-10:.01:2]; cbarxticks = [-30:.5:10];
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');

%%

if 1==2 % if add nonlinear drag contours
figure(100)
hold on;
cax    = caxis;
uo     = 0;             % amplitude of mean nontidal vel
utidal = (2*kineticEnergyAv_TOT).^1/2; % amplitude of mean tidal vel
tilomhere = 1;
varp2  = (abs(tilomhere.*tilcesq_TOT)) ./ ( tilT_TOT.*(uo^2 + utidal).^(1/2) );
%contour(tilT_TOT,tilcesq_TOT,log10(varp2),[1:.005:5],'k:');
cv = contour(tilT_TOT,tilcesq_TOT,log10(varp2),20,'k:'); %clabel(cv)
caxis(cax);
hold off
print([workdir,'/Figs/',figvar,'_',component,compsub,'_wCdcont'],'-dpng');
end


%%

figvar    = 'energy'; % energy density (not multiplied by thickness)
varp      =  kineticEnergyAv_TOT + sign(tilcesq_TOT).*potentialEnergyAv_TOT;
varpnayme = 'log$_{10}$(total energy density)';
cl = [-10:.1:4]; cbarxticks = [-15:1:15];
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');


%%

figvar    = 'kineticEnergy'; % kinetic energy (not multiplied by thickness)
varp      = kineticEnergyAv_TOT; 
varpnayme = 'log$_{10}$(kinetic energy density)';
cl = [-10:.1:4]; cbarxticks = [-15:1:15];
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');


%%


figvar    = 'potentialEnergy'; % potential energy (not multiplied by thickness)
varp      = sign(tilcesq_TOT).*potentialEnergyAv_TOT;
varpnayme = 'log$_{10}$(potential energy density)';
cl = [-10:.1:4]; cbarxticks = [-15:1:15];
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');
  

%%


figvar    = 'abs(knFsF)'
varp      = abs(knFsF_TOT); 
varpnayme = 'abs(knFsF)';
cl = [0:.01:1]*1.5; 
cbarxticks = [-1:.1:1]*1.5;
%
figure(100);clf
contourf(tilT_TOT,tilcesq_TOT,(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');
 
%%

figvar    = 'phaseknFsF'
varp      = angle(knFsF_TOT)/(2*pi); 
varpnayme = 'phase(knFsF)';
if min(varp(:)) < 0 & max(varp(:)) <= 0,          % case all lags are negative
    cl =[-0.5:.001:0];   cbarxticks = [-1:.05:1];
elseif min(varp(:)) < 0 & max(varp(:)) > 0,       % special case of G20 for which response can have pos and neg lags
    cl =[-0.5:.001:0.5];  cbarxticks = [-1:.05:1];
else 
    cl =[0:.001:0.5];    cbarxticks = [-1:.05:1]; % case all lags are positive
end

figure(200);clf
contourf(tilT_TOT,tilcesq_TOT,(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
cbar = colorbar('Eastoutside'); cbar.Label.Interpreter = 'latex'; cbar.Label.String = varpnayme;
xlabel('${\tilde T}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log'); set(gca,'xscale','log');  
set(gca,'fontsize',fsize); set(gca,'Xtick',10.^[-10:1:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component,compsub],'-dpng');




%%


end % end of loop for negtilcesq_sw 

