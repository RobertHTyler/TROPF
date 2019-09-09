function field = SH2Field(Cns, lon,lat,t, s,om,phase0, swCC)
% function field = SH2Field(Cns, lon,lat,t, s,om,phase0, swCC)
% 
% Map Cns SH coefs to lon,lat grid to produce field at time t.
% 
% In:
% Cns:              column vector of SH coefs
% lon,lat:          2D arrays of longitude and of latitude (all grid points; i.e. number of elements in array = number of grid points)
% t:                time (scalar)
% s:                order/rank (scalar)
% om:               'omega', the temporal frequency (a scalar which may be positive or negative)
% phase0:           arbitrary phase offset to add to phase argument
%
% Out: 
% field:            2D array (lon,lat) of mapped Cns coefs
%
% Notes:
% 1) Assumes lon and lat are matrices indexed with lon running first dim and
% lat the second (i.e. [Nlon,Nlat] = size(lon)) 
% 2) om (omega) and t (time) can be dimensional or nondimensional as they only appear in this functions as the product om*t.
%
% R. Tyler (5 Jan., 2019)
%
%%

N = length(Cns); % number of SH terms in expansion

% Create array of spherical-harmonic functions:
var_Ynsom_arr = Ynsom_arr(lon,lat,t, N, s,om,phase0);

% Weight by coefs Cns and reshape: 
[Nlon,Nlat] = size(lon); % grid dimensions
field       = permute( reshape( var_Ynsom_arr*Cns , Nlat,Nlon ), [2 1] );
field       = field + swCC*conj(field); 

