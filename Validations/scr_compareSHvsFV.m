% Validation script for comparing tropf.m results with those from a
% finite-volume model
%
% R. Tyler, 10 March, 2019
% %
%% % % % Housekeeping:

clear
TROPFPATH = '~/Desktop/TROPF';        % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;     % run startup script to set paths to TROPF libraries
workdir = [TROPFPATH,'/Validations']; % working directory
cd(workdir);                          % move to working directory




%% % % % FV METHOD % % % %:

% % % Define grid:
DELLON = 1; DELLAT = DELLON; 
[lon,lat]   = ndgrid(DELLON:DELLON:DELLON+(360-DELLON), -90:DELLAT:90); 
[Nlon,Nlat] = size(lon); Ngrid = Nlon*Nlat; % grid dimensions
colat       = 90 - lat;     % colatitude (deg)


% % % Define input parameters:

% edit these physical and configuration parameters:
n=2; s=2; omt = 0;              % degree, order n,s of forcing; omt = temporal phase (omega*t)
tilG   =   (legendrens(n,s,cos(colat*pi/180)).*(exp(1i*(s*lon*pi/180 -omt)))); % forcing potential tilG
tilGscalefac = max(abs(tilG(:))); tilG = tilG/tilGscalefac;   % normalize for amplitude=1
tilOm   =  1; 
tilf    =  tilOm*sin(lat*pi/180);  % nondim. rotation rate
tilom   = -0.766; 
PT      = -1i*tilom;
tilalp  =  1e-1;
tilcesq =  1e-2;

% operators:
Ltilalp =  tilalp;
Lh      =  1;
LV      =  1./tilcesq;
 
% % define PDE parameters: 
Dcoef  =  Lh.*(PT+tilalp)./( (PT+tilalp).^2 + tilf.^2 +eps); % nondim. complex diffusion coefficient
%Dcoef  =  Lh./( (PT+tilalp) + tilf.^2./(PT+tilalp) );   % nondim. complex diffusion coefficient
Psis   =  Lh.*tilf./( (PT+tilalp).^2 + tilf.^2 ) ;       % nondim. complex string function
[Dxvar,Dyvar,Gradstruct] = gradGlobe(lon,lat,1,  Psis ); % nondim. gradient of string function
Vx   = Dyvar; Vy = -Dxvar;                               % nondim. complex string velocity
Hcoef  =  PT*LV;                                         % nondim. helmholtz coef.
forc   = -PT*LV.*tilG;                                   % nondim. forcing term


% % % Build coef matrix and solve:
[metrik]      =  operatorFVadvdiff(lon,lat,1,Dcoef,Hcoef,Vx,Vy);
tilpMtilG_FV  = (reshape( metrik.L \ (forc(:).*metrik.vol(:)),  Nlon,Nlat));
%tilpMtilG_FV = (1/2)*( reshape( metrik.L \ (forc(:).*metrik.vol(:)),  Nlon,Nlat)...
%              + reshape( conj(metrik.L) \ conj(forc(:).*metrik.vol(:)),  Nlon,Nlat) );
tilp_FV       = tilpMtilG_FV + tilG;% (1/2)*(tilG+conj(tilG));



%tilpMtilG_FV = cos(lat*pi/180);
% velocities
r=1;
tilFpx = 0; tilFpy = 0; 
[junkx,junky,Gradstruct] = gradGlobe(lon,lat,r,  tilpMtilG_FV );
u = (PT + tilalp)./((PT+tilalp).^2+tilf.^2) .* (-junkx + tilFpx)...
  + ( tilf./((PT+tilalp).^2+tilf.^2)) .* (-junky + tilFpy); 
v = (PT + tilalp)./((PT+tilalp).^2+tilf.^2) .* (-junky + tilFpy)...
  + (-tilf./((PT+tilalp).^2+tilf.^2)) .* (-junkx + tilFpx); 

%u = junkx; v = junky;
u(:,[1,Nlat]) = 0; v(:,[1,Nlat]) = 0;
u([1 Nlon],:)=0; v([1 Nlon],:)=0;

u = real(u); v = real(v);

figure(50); clf;% subplot(211)
%contour(lon,lat,tilp_FV,20); colorbar
%contour(lon,lat,tilpMtilG_FV,20); colorbar
contour(lon,lat,real(tilG),20); colorbar
hold on; quiver(lon,lat,u,v,'r'); hold off
%hold on; quiver(lon,lat,udE+urE,udN+urN,'k'); hold off
axis tight



%% % % % SPHERICAL HARMONIC METHOD % % % %:
tilt    = 0; % time zero 
phase0  = 0; % 
N       = 500;
tilalpd = tilalp(1); 
tilalpr = tilalpd; tilalpp = 0;
Gns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1); Kns = zeros(N,1);
Gns(n-s+1) = 1i*(1/2) * max(abs(legendrens(n,s,[-1:1e-3:1]))).^(-1); % amplitude of force for non-dim calcs (for |mathfrakG|=1 m)
tilnusqns  =  (1 + 1i*tilalpp/tilom)./tilcesq;

% solve:
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);

Gfield = SH2Field(-1i*Gns,     lon,lat,tilt, s,tilom,phase0,1);
Dfield = SH2Field( Dns,        lon,lat,tilt, s,tilom,phase0,1);
Rfield = SH2Field(-1i*Rns,     lon,lat,tilt, s,tilom,phase0,1);
pfield = SH2Field(-1i*pns,     lon,lat,tilt, s,tilom,phase0,1);
%
udE =  SH2Field_gradE( Dns,    lon,lat,tilt, s,tilom,phase0,1);
udN =  SH2Field_gradN( Dns,    lon,lat,tilt, s,tilom,phase0,1);
urE =  SH2Field_gradN(-1i*Rns, lon,lat,tilt, s,tilom,phase0,1);
urN = -SH2Field_gradE(-1i*Rns, lon,lat,tilt, s,tilom,phase0,1);
%
uE = udE + urE; uN = udN + urN; 
uE(:,[1,end]) = 0; uN(:,[1,end]) = 0; 

PTpfield  = SH2Field((-1i*tilom)*(-1i*pns),     lon,lat,tilt, s,tilom,phase0,1);
workfield = (1/tilcesq)*Gfield.*PTpfield;

EKfield      = 0.5*(uE.^2 + uN.^2);
dissipfield  = 2*tilalpd*EKfield;

%%

figure(30); clf; subplot(311)
contour(lon,lat,real(tilp_FV),50); colorbar
%hold on; quiver(lon,lat,u,v); hold off
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Finite Volume');
axis tight

figure(30); subplot(312)
contour(lon,lat,pfield,50); colorbar
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Spherical Harmonic');
axis tight
%hold on; quiver(lon,lat,uE,uN); hold off
axis tight

figure(30); subplot(313)
contour(lon,lat,real(pfield-tilp_FV),50); colorbar
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Spherical Harmonic minus Finite Volume');
%hold on; quiver(lon,lat,urE,urN); hold off
axis tight
colormapBlueToRed
orient tall; print -dpng ./valFVvsSH; orient portrait

mabs(tilp_FV,pfield)


% 
% %%
% velamp = sqrt( (real(udE+urE)).^2 +(real(udN+udN)).^2 );
% 
% figure(50); clf;% subplot(211)
% contour(lon,lat,velamp,20); colorbar
% hold on; quiver(lon,lat,u,v,'r'); hold off
% hold on; quiver(lon,lat,udE+urE,udN+urN,'k'); hold off
% axis tight
% colormapBlueToRed
% print -dpng valSHvsFV
% 
% uE = udE+urE; uN = udN+urN;
% 
% mabs( u, udE+urE ), mabs( v, udN+urN ),
% 
% %%
% 
% figure(31); clf; %subplot(211)
% contour(lon,lat,Dfield,50); colorbar
% hold on; quiver(lon,lat,udE,udN); hold off
% hold on; quiver(lon,lat,u,v,'m'); hold off
% axis tight
% %%
% if 1==2
% figure(31); subplot(212)
% contour(lon,lat,Rfield,50); colorbar
% hold on; quiver(lon,lat,urE,urN); hold off
% axis tight
% end
%  