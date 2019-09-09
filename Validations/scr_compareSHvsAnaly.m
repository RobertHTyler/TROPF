% Validation script to compare tropf.m vs Analytical results

% R. Tyler 1 March 2019


%% % % Housekeeping:

clear
TROPFPATH = '~/Desktop/TROPF';       % set path of TROPF root directory
cd(TROPFPATH); scr_startup_tropf;    % run startup script to set paths to TROPF libraries
workdir = [TROPFPATH,'/Validations'];   % working directory
cd(workdir);                         % move to working directory



%%
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
tilOm   =  0; 
tilf    =  tilOm*sin(lat*pi/180);  % nondim. rotation rate
tilom   = -0.766; 
%tilom   = -0.766*0 + eps; 
PT      = -1i*tilom;
%tilalp  =  1e-1;
tilalp  =  0;
tilcesq =  1e-2;

% operators:
Ltilalp =  tilalp;
Lh      =  1;
LV      =  1./tilcesq;
 

tilt    = 0; % time zero 
phase0  = 0; % 
N       = 500;
tilalpd = tilalp(1); 
tilalpr = tilalpd; tilalpp = 0;
Gns = zeros(N,1); dns = zeros(N,1); ens = zeros(N,1); Kns = zeros(N,1);
Gns(n-s+1) = 1i*(1/2) * max(abs(legendrens(n,s,[-1:1e-3:1]))).^(-1); % amplitude of force for non-dim calcs (for |mathfrakG|=1 m)
tilnusqns  =  (1 + 1i*tilalpp/tilom)./tilcesq;



%% % Solve using tropf.m:
[Dns,Rns,pns, calWns,calDns,calEKns,calEPns,knFsF] = tropf(N,tilOm, tilom,s,Gns,Kns,dns,ens, tilalpd,tilalpr,tilnusqns);
pfield = SH2Field(-1i*pns,     lon,lat,tilt, s,tilom,phase0,1);
Gfield = SH2Field(-1i*Gns,     lon,lat,tilt, s,tilom,phase0,1);
 

%% % Analytical solution:

nF = n; 
%tilCf = (-1i*tilom + Ltilalp) ./ ( (-1i*tilom + Ltilalp)^2 );
%tilCf = 1 ./ ( (-1i*tilom + Ltilalp) );

%tilCf  =  Lh*(PT+tilalp)./( (PT+tilalp).^2 + tilf.^2); % nondim. complex diffusion coefficient


tilCf  =  1./(PT+tilalp); % nondim. complex diffusion coefficient

%tilp_analy =  ( 1./(tilCf*nF*(nF+1)/(1i*tilom*LV) -1)  + 1 ) .* Gfield;  
%tilp_analy = ( (n*(n+1)*tilcesq/tilom^2 -1)^(-1) + 1) * Gfield;



tilp_analy = real(  ( (tilCf*n*(n+1)*tilcesq/(1i*tilom) -1)^(-1) + 1) * Gfield  );

 
%tilp_analy = ( (tilCf*(-1i*tilom))*(n*(n+1)*tilcesq/tilom^2 -1)^(-1) + 1) * Gfield;



figure(30); clf; subplot(311)
contour(lon,lat,real(tilp_analy),50); colorbar
%hold on; quiver(lon,lat,u,v); hold off
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Analytical');
axis tight

figure(30); subplot(312)
contour(lon,lat,pfield,50); colorbar
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Spherical Harmonic');
axis tight
%hold on; quiver(lon,lat,uE,uN); hold off
axis tight

figure(30); subplot(313)
contour(lon,lat,(pfield - tilp_analy),50); colorbar
xlabel('longitude (degrees)'); ylabel('latitude (degrees)'); title('Spherical Harmonic minus Analytical');
%hold on; quiver(lon,lat,urE,urN); hold off
axis tight
colormapBlueToRed


orient tall; print -dpng ./valAnalyvsSH_0tilalp; orient portrait
%orient tall; print -dpng ./valAnalyvsSH_0tilom; orient portrait

mabs(tilp_analy,pfield)


% Note maximum residuals for tilalp=0 case:
% 3.7470e-16 (N=500)
% 1.6653e-16 (N=1000)


