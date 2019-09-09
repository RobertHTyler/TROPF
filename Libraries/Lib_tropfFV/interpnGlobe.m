function ZOut = interpnGlobe(lonIn,latIn,ZIn, lonOut,latOut,method)
%
% Interpolation of Zin over continuous 2D spherical surface with
% coordinates lonIn,latIn (degrees). This scheme wraps values across
% longitudinal ends, and extrapolates poleward using longitudinal average
% of most poleward data.
%
% Assumes:
% 1) regular grid: first dim is monontonically increasing longitude, second
% is monotonically increasing latitude; no overlap points across wrap longitude
% (0/360); (assumes  0< lonIn <= 360)
% 2) Does not (yet) weight for decreasing surface area of poleward pixel
% areas (i.e. proximity of neighbors is by angular degrees rather than meters.) 
% 3) for volumetric weighting, assumes r = radius of Earth
%
% Default method = 'spline'; 

if nargin == 5, method = 'linear'; end
    
    

volumetric_weighting = 0; 

%%% Metriks of input params:
Rearth=6371e3; r=Rearth;
lon=lonIn; lat=latIn; 
[Mm,Nn]=size(lon);   % size of array
% Determine resolution:
% assumes resolution of center of domain and (so far) that grid is uniform!
DELi=(lon(round(Mm/2),1)-lon(round(Mm/2)-1,1)); 
DELj=(lat(1,round(Nn/2))-lat(1,round(Nn/2)-1));                                                           
% Define some preliminary parameters:						    
sinlat=sin(lat*pi/180);
sinlatph=sin((lat+DELj/2)*pi/180);
sinlatmh=sin((lat-DELj/2)*pi/180);
coslat=cos(lat*pi/180);
coslatph=cos((lat+DELj/2)*pi/180);
coslatmh=cos((lat-DELj/2)*pi/180);
delx=r*coslat*DELi*pi/180;
dely=r.*DELj*pi/180+0*lon;       
% Define surf "areas" of sides of each cell, and "vol" of cell:
Aw=DELj*pi/180*r+0*lon;
Ae=Aw;
As=DELi*pi/180*coslatmh*r;
An=DELi*pi/180*coslatph*r;
vol=DELi*pi/180 * (sinlatph-sinlatmh).*r.^2; %volume of cell (i.e. surface area ftc)
%vol=delx.*dely; % 
volIn = vol; % this is the volume (area in 2D) param we need for volumetric averaging


if volumetric_weighting==1
    ZIn = ZIn.*volIn;
end

% Extend matrices for wrap in longitude (assumes  0<= lonIn <= 360): 
nlap = 10; 
lonIn_extend = [lonIn(end-nlap:end,:)-360; lonIn; lonIn(1:nlap,:)+360];
latIn_extend = [latIn(end-nlap:end,:); latIn; latIn(1:nlap,:)];
ZIn_extend = [ZIn(end-nlap:end,:); ZIn; ZIn(1:nlap,:)];
vol_extend = [vol(end-nlap:end,:); vol; vol(1:nlap,:)];
% further extend to extend input data to South Pole: 
minlat = min(latIn(:)); 
numofmissinglatsS = ceil((minlat-(-90))/DELj);
for ii=1:numofmissinglatsS
lonIn_extend = [lonIn_extend(:,1), lonIn_extend];
latIn_extend = [latIn_extend(:,1)-DELj, latIn_extend];
ZIn_extend   = [0*ZIn_extend(:,numofmissinglatsS+1)+mean(ZIn_extend(:,numofmissinglatsS+1)), ZIn_extend];
end
% further extend to extend input data to North Pole: 
maxlat = max(latIn(:)); 
numofmissinglatsN = ceil((90-maxlat)/DELj);
for ii=1:numofmissinglatsN
lonIn_extend = [lonIn_extend, lonIn_extend(:,end)];
latIn_extend = [latIn_extend, latIn_extend(:,end)+DELj];
ZIn_extend   = [ZIn_extend, 0*ZIn_extend(:,end-numofmissinglatsN+1)+mean(ZIn_extend(:,end-numofmissinglatsN+1))];
end

% Point interp:

ZOut = interpn(lonIn_extend,latIn_extend, ZIn_extend, lonOut,latOut,method); 
%indd = ~isnan(ZIn_extend);
%ZOut = griddata(lonIn_extend(indd),latIn_extend(indd), ZIn_extend(indd), lonOut,latOut,'linear'); 






if volumetric_weighting == 1

%% metriks of Output params:
Rearth=6371e3; r=Rearth;
lon=lonOut; lat=latOut; 
[Mm,Nn]=size(lon);   % size of array

% Determine resolution:
% assumes resolution of center of domain and (so far) that grid is uniform!
DELi=(lon(round(Mm/2),1)-lon(round(Mm/2)-1,1)); 
DELj=(lat(1,round(Nn/2))-lat(1,round(Nn/2)-1)); 
                                                           
% Define some preliminary parameters:						    
sinlat=sin(lat*pi/180);
sinlatph=sin((lat+DELj/2)*pi/180);
sinlatmh=sin((lat-DELj/2)*pi/180);
coslat=cos(lat*pi/180);
coslatph=cos((lat+DELj/2)*pi/180);
coslatmh=cos((lat-DELj/2)*pi/180);
delx=r*coslat*DELi*pi/180;
dely=r.*DELj*pi/180+0*lon;
       
% Define surf "areas" of sides of each cell, and "vol" of cell:
Aw=DELj*pi/180*r+0*lon;
Ae=Aw;
As=DELi*pi/180*coslatmh*r;
An=DELi*pi/180*coslatph*r;
volOut=DELi*pi/180 * (sinlatph-sinlatmh).*r.^2; %volume of cell (i.e. surface area ftc)
%vol=delx.*dely; % 



ZOut=ZOut./(volOut+eps) ; 

ZOut(:,1)=mean(ZOut(:,2)); ZOut(:,end)=ZOut(:,end-1); 
end


%% Check with plot:
if 1==2
addpath('/Users/rtyler/Desktop/Tidal_mag_forward/lib_Tyler/m_map')  
msz=5; % marker size

figure(1002); clf; cla; 
m_proj('stereographic','lat',[-90],'long',0,'radius',[45]);
m_contourf(lonOut,latOut,ZOut,30); shading flat;  axis off
m_coast('color','m','linestyle','.','markersize',msz)


figure(1000); clf; cla; 
m_proj('mercator')%,'lat',[-90],'long',0,'radius',[45]);
m_contourf(lonOut,latOut,ZOut,180); shading flat;  axis off
m_coast('color','m','linestyle','.','markersize',msz)

figure(1001); clf; cla; 
m_proj('stereographic','lat',[90],'long',180,'radius',[45]);
m_contourf(lonOut,latOut,ZOut); shading flat;  axis off
m_coast('color','m','linestyle','.','markersize',msz)
end

