% Plot various solution parameters in the solution domain. This assumes
% that the solutions have been built (e.g. using the  ./buildParamSpace_*.m file) and stored in the ./Results directory.
% 


clear
TROPFPATH  = '~/Desktop/TROPF';      % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;   % run startup script to set paths to TROPF libraries
resultsdir = [TROPFPATH,'/Studies/NonSynchRotation/ArbRotation/Results'];  % working directory
workdir    = [TROPFPATH,'/Studies/NonSynchRotation/ArbRotation'];  % working directory
cd(workdir);                        % move to working directory

fsize = 16;

 
%% Choose solution set to load

component = 'G22W_RD';
% component = 'G22WN_RD';
% component = 'G21W_RD';
% component = 'G21WN_RD';
%  
% component = 'G22W_NR';
% component = 'G22WN_NR';
% component = 'G21W_NR';
% component = 'G21WN_NR';
 
load([resultsdir,'/',component]);

%% Check some equivalences:

mabs(workAv_TOT,dissipationAv_TOT) % should be small since average work and dissipation are equal

%mabs(dissipationAv_TOT,  2*tilalpd*kineticEnergyAv_TOT)


% for tilom near zero, a roundoff imaginary part might be added because the
% real part is so small. Check 
max(abs( imag(workAv_TOT(:))  )), % this should be small
workAv_TOT = real(workAv_TOT); 

%% Plots:


%%
figvar    = 'power'; % power density 
varp      = workAv_TOT; 
varpnayme = 'log$_{10}({\tilde { {\cal P }} }) $';
cl = [-10:.05:2]; cbarxticks = [-10:.5:10];
%
figure(100);clf
contourf(tilom_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde \omega}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex';
cbar.Label.String = varpnayme;
xlabel('${\tilde \omega}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log');   
set(gca,'fontsize',fsize); set(gca,'Xtick',[-10:0.5:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component],'-dpng');

 

%%
figvar    = 'tilcesqpower'; % power density multiplied by tilcesq (prop. to layer thickness)
varp      = abs(tilcesq_TOT).*workAv_TOT; 
varpnayme = 'log$_{10}(|{\tilde c_e}^2|{\tilde { {\cal P }} }) $';
cl = [-10:.05:2]; cbarxticks = [-10:.5:10];
%
figure(100);clf
contourf(tilom_TOT,tilcesq_TOT,log10(varp),cl,'linestyle','none'); shading flat; 
xlabel('${\tilde \omega}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
cbar = colorbar('Eastoutside','xtick',cbarxticks); cbar.Label.Interpreter = 'latex';
cbar.Label.String = varpnayme;
xlabel('${\tilde \omega}$','Interpreter','latex'); ylabel('${\tilde c_e}^2$','Interpreter','latex'); 
ylabel('${\tilde c_e}^2$','Interpreter','latex'); set(gca,'yscale','log');   
set(gca,'fontsize',fsize); set(gca,'Xtick',[-10:0.5:10])
if tilcesq_TOT(1) > 0, set(gca,'Ytick',10.^[-10:1:10]); else, set(gca,'Ytick',fliplr(-10.^[-10:1:10])); end; 
colormap jet; caxis([min(cl) max(cl)]); grid
print([workdir,'/Figs/',figvar,'_',component],'-dpng');

 