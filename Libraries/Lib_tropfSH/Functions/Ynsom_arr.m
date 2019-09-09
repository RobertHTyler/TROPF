function Ynsom_arr = Ynsom_arr(lon,lat,t, N, s,om,phase0)
% function Ynsom_arr = Ynsom_arr(lon,lat,t, N, s,om,phase0)
%
% Returns the 2D array (matrix) where each row is the range of spherical
% harmonics over degree at a specific lon,lat,t. The columns span the various lon,lat pairs. 
%
% In:
% lon,lat:          2D arrays of longitude/latitude pairs
% t:                time (scalar)
% N:                Number of elements in spherical-harmonic expansion
% s:                order/rank (scalar)
% phase0:           arbitrary phase offset to add to phase argument
% om:               'omega', the temporal frequency (a scalar which may be positive or negative)
%
% Notes:
% Assumes lon and lat are matrices indexed with lon running first dim and
% lat the second (i.e. [Nlon,Nlat] = size(lon)) 
%
% t and om may be dimensional or nondimensional (i.e. tilt, tilom); they only appear as the
% product om*t here.)
%
% R. Tyler (5 Jan., 2019)
%%

[Nlon,Nlat] = size(lon); Ngrid = Nlon*Nlat; % grid dimensions
colat       = 90 - lat                    ; % colatitude (deg)
x           = permute(cos(pi/180*colat(1,:)),[2 1]); % col. vec. of cos(colat) 
Ntrunc      = N + s - 1                   ; % truncation degree of SH expansion
nvec        = permute([s:Ntrunc],[2 1])   ; % col. vec. of SH degrees involved

% Create Associated Legendre function at each lat, then replicate through lons to produce
% full array:
Pnsx =  zeros(Nlat,length(nvec));
for ii = 1:length(nvec);
    Pnsx(:,ii) = Pns_xvec(nvec(ii),s,x);
end
Pnsx_arr = repmat(Pnsx,Nlon,1); % Ngrid X N matrix ( first column corresponds to locations of all lats and first lon, followed by all lats and second lon etc.)

% Create phase at all locations and then replicate through degrees to produce full
% array:
lons_vec   = permute(lon,[2 1]); lons_vec = lons_vec(:);
phayse_vec = exp(1i*(s*(pi/180*lons_vec) - om*t + phase0));       % phase at all locations
phayse_arr = repmat(phayse_vec,1,length(nvec));

% Full array:
Ynsom_arr  = Pnsx_arr.*phayse_arr;


