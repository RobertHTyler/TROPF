function [metrik] = operatorFVadvdiff(lon,lat,r,Dcoef,Hcoef,Vx,Vy);
% function [metrik] = operatorFVadvdiff(lon,lat,r,Dcoef,Hcoef,Vx,Vy);
%
% This function returns the finite-volume matrix operator and metrics for
% the cell-volume integrated advection-diffusion equation:
%
% integral(   divergence( Dcoef(gradient(var)) + var*Uvec )  + Hcoef*var  ) = integral(  force  )
%
% where var is the variable to be solved for and force is the prescribed
% force (these variables do not appear in this function because only the
% operator is being built here). The input and outputs are as follows:
%
% INPUT:
%
% lon, lat      : longitude, latitude; (degrees); 2-D arrays over complete sphere; no repeating longitude values; Can go to +- 90 degrees latitude
% r             : radius (m);
% Dcoef         : Diffusion coefficient; 
% Hcoef         : coefficient
% Vx,Vy         : Advection velocities in 'eastward' and 'northward' directions
%
%
% OUTPUT:
%
% metrik.Aw      : area/length of west wall of cell
% metrik.Ae      : area/length of east wall of cell
% metrik.As      : area/length of south wall of cell
% metrik.An      : area/length of north wall of cell
% metrik.Vol     : volume/surface area of cell
% metrik.Lfluxes : outward fluxes through cell walls
% metrik.L       : final operator associated with equation above




% R. Tyler 10 Dec. 2018

%%
display('starting opFVadvdiff')

%%% Some renaming of input params:

% If Dcoef is entered as a structure, then expand:
if isstruct(Dcoef),
Dcoefw=Dcoef.w; Dcoefe=Dcoef.e; Dcoefs=Dcoef.s; Dcoefn=Dcoef.n; 
else % else create Dcoef at the four cell walls by averaging:
Dcoefw=avhw(Dcoef); Dcoefe=avhe(Dcoef); Dcoefs=avhs(Dcoef); Dcoefn=avhn(Dcoef);
Dcoefw(1,:)=0.5*(Dcoef(1,:)+Dcoef(end,:)); Dcoefe(end,:)=Dcoef(1,:); % adjust for wrap points
end
% If Vx is entered as a structure, then expand:
if isstruct(Vx),
Vxw=Vx.w; Vxe=Vx.e; 
else % else create Vx at the W/E cell walls by averaging:
Vxw=avhw(Vx); Vxe=avhe(Vx);
Vxw(1,:)=0.5*(Vx(1,:)+Vx(end,:)); Vxe(end,:)=Vxw(1,:);% adjust for wrap points
end
% If Vy is entered as a structure, then expand:
if isstruct(Vy),
Vys=Vy.s; Vyn=Vy.n; 
else % else create Vy at the S/N cell walls by averaging:
Vys=avhs(Vy); Vyn=avhn(Vy); 
end


% Determine resolution (assumes grid is uniform!)
[Nlon,Nlat] = size(lon); 
Ngrid = Nlon*Nlat;
DELi = lon(2,1)-lon(1,1); DELj = lat(1,2)-lat(1,1);

% Find boundary indices:
bw = find(lon==min(lon(:))); % west boundary
be = find(lon==max(lon(:))); % east boundary
bs = find(lat==min(lat(:))); % south boundary
bn = find(lat==max(lat(:))); % north boundary
 


                                                          
% Define some preliminary parameters:	     
sinlat   = sin(lat*pi/180);
sinlatph = sin((lat+DELj/2)*pi/180);
sinlatmh = sin((lat-DELj/2)*pi/180);
coslat   = cos(lat*pi/180);
coslatph = cos((lat+DELj/2)*pi/180);
coslatmh = cos((lat-DELj/2)*pi/180);
delx     = r*coslat*DELi*pi/180;
dely     = r.*DELj*pi/180 + zeros(Nlon,Nlat);
       
% Define surf "areas" of sides of each cell, and "vol" of cell:
Aw  = DELj*pi/180*r + zeros(Nlon,Nlat);
Ae  = Aw;
As  = DELi*pi/180*coslatmh*r;
An  = DELi*pi/180*coslatph*r;
vol = DELi*pi/180*(sinlatph-sinlatmh).*r.^2; %volume of cell (i.e. surface area for 2D case)


% % Build matrices: 

% Build diffusive fluxes out of cell:
Cw = -Dcoefw.*Aw./avhw(delx); Cw = Cw(:);
Ce = -Dcoefe.*Ae./avhe(delx); Ce = Ce(:);
Cs = -Dcoefs.*As./avhs(dely); Cs = Cs(:);
Cn = -Dcoefn.*An./avhn(dely); Cn = Cn(:);
difFluxw = spdiags([-Cw(:),[Cw(2:end);nan]],                  [0, -1],    Ngrid,Ngrid);
difFluxe = spdiags([-Ce(:),[nan;Ce(1:end-1)]],                [0, 1],     Ngrid,Ngrid);
difFluxs = spdiags([-Cs(:),[Cs(1+Nlon:end);nan*ones(Nlon,1)]],[0, -Nlon], Ngrid,Ngrid);
difFluxn = spdiags([-Cn(:),[nan*ones(Nlon,1);Cn(1:end-Nlon)]],[0, Nlon],  Ngrid,Ngrid);
 
% Build advective fluxes out of cell:
aCw = Vxw.*Aw /2; 
aCe = Vxe.*Ae /2; 
aCs = Vys.*As /2; 
aCn = Vyn.*An /2; 
aCw(1,:)   = 0.5*(Vxw(1,:)+Vxw(end,:)) .*Aw(1,:) *(1/2); 
aCe(end,:) = aCe(1,:); 
aCw = aCw(:); aCe = aCe(:); aCs = aCs(:); aCn = aCn(:);
advFluxw = spdiags([-aCw(:),-[aCw(2:end);nan]],                  [0, -1],    Ngrid,Ngrid);
advFluxe = spdiags([ aCe(:), [nan;aCe(1:end-1)]],                [0, 1],     Ngrid,Ngrid);
advFluxs = spdiags([-aCs(:),-[aCs(1+Nlon:end);nan*ones(Nlon,1)]],[0, -Nlon], Ngrid,Ngrid);
advFluxn = spdiags([ aCn(:), [nan*ones(Nlon,1);aCn(1:end-Nlon)]],[0, Nlon],  Ngrid,Ngrid);


% Total Fluxes:
Fluxw = difFluxw + advFluxw; 
Fluxe = difFluxe + advFluxe; 
Fluxs = difFluxs + advFluxs; 
Fluxn = difFluxn + advFluxn; 



% Build other terms:
integHcoef = spdiags(Hcoef(:).*vol(:), 0, Ngrid,Ngrid);



% 
% % % Adjust for polar boundary:
% Fluxs(bs,:) = 0; Fluxn(bn,:) = 0;
% 
% Cjunk=zeros(Nlon,Nlat)(:) + full(mean(abs(diag(Fluxw))));
% 
% Cjunk=1e6;
% 
% if 1==2
%  % north pole:
% bjunk=find(lon>min(lon(:)) & lat==90); 
%  Fluxn=Fluxn+sparse(bjunk,bjunk,-Cjunk(bjunk),Ngrid,Ngrid);
%  Fluxn=Fluxn+sparse(bjunk,bjunk-1,Cjunk(bjunk),Ngrid,Ngrid);
% % set one polar point equal to average of equatorward neighbors:
% bjunk=find(lon==min(lon(:)) & lat==90); 
%  Fluxn=Fluxn+sparse(bjunk,bjunk,-Cjunk(bjunk),Ngrid,Ngrid);
%  for ijunk=1:Nlon
%  Fluxn=Fluxn+sparse(bjunk,bjunk+(-Nlon-1+ijunk),Cjunk(bjunk)/Nlon,Ngrid,Ngrid);
%  end
% end
% 
% 
% if 1==1 % this BC instead should set value to zero at SP:
%  % south pole:
% bjunk=find(lat==-90); 
% Fluxs(bjunk,bjunk) = ones(length(bjunk))
% end
%  
% 
% if 1==1 % this BC instead should set value to zero at NP:
%  % north pole:
% bjunk=find(lat==90); 
%  Fluxs=Fluxs+sparse(bjunk,bjunk,-Cjunk(bjunk),Qq,Qq);
% end
% 



% % Adjust for polar boundary:
%Fluxs(bs,:) = 0; Fluxn(bn,:) = 0;

%  % south pole:
%bjunk=find(lat==-90); integHcoef(bjunk,bjunk) = 1; 

%  % north pole:
%bjunk=find(lat==90); integHcoef(bjunk,bjunk) = 1; 


% Adjust flux for E/W wrap:
% flux through west side of cell:
bjunk=find(lon==min(lon(:)) & lat>min(lat(:)) & lat<max(lat(:))); 
Fluxw(bjunk,:) = 0;
Fluxw = Fluxw + sparse(bjunk,bjunk,        -Cw(bjunk), Ngrid,Ngrid);
Fluxw = Fluxw + sparse(bjunk,bjunk-1+Nlon,  Cw(bjunk), Ngrid,Ngrid);
% flux through east side of cell:
bjunk=find(lon==max(lon(:)) & lat > min(lat(:)) & lat<max(lat(:)));
Fluxe(bjunk,:) = 0;
Fluxe = Fluxe + sparse(bjunk,bjunk,       -Ce(bjunk), Ngrid,Ngrid);
Fluxe = Fluxe + sparse(bjunk,bjunk+1-Nlon, Ce(bjunk),Ngrid,Ngrid);



metrik.Aw = Aw; metrik.Ae = Ae; metrik.As = As; metrik.An = An; 
metrik.vol = vol; 

metrik.Lfluxes =  Fluxw + Fluxe + Fluxs + Fluxn ;
metrik.L       =  metrik.Lfluxes + integHcoef;





