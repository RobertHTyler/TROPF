function [fieldMov, fieldMov_gradE,fieldMov_gradN] = buildMovie(Cns, lon,lat,tiltVec, s,tilom,phase0, swCC)
% function [fieldMov, fieldMov_gradE,fieldMov_gradN] = buildMovie(Cns,
% lon,lat,tiltVec, s,tilom,phase0, swCC) 
%
% This function projects spherical harmonic coefs (Cns) onto lon,lat grid
% at times tiltVec
%
% IN:
% Cns:              column vector of SH coefs to map in 3D (lon,lat,tymes)
% lon,lat:          2D arrays of longitude and of latitude (all grid points; i.e. number of elements in array = number of grid points)
% tiltVec:          time (vector); times for movie frames
% s:                order/rank (scalar)
% tilom:            'tilde omega', the temporal frequency (a scalar which may be positive or negative)
% phase0:           arbitrary phase offset to add to phase argument
% swCC:             swCC = 1, unless the field expanded has a 1i* multiplier, in which case swCC = -1 (because we need to subtract rather than add the complex conjugate in the latter case)
%
% OUT: 
% fieldMov:         3D array (lon,lat,tiltVec) of mapped Cns coefs
% field_gradE:      3D array (lon,lat,tiltVec) of mapped 'eastward'  gradient component of fieldMov
% field_gradN:      3D array (lon,lat,tiltVec) of mapped 'northward' gradient component of fieldMov
%
% Notes:
% 1) Assumes lon and lat are matrices indexed with lon running first dim and
% lat the second (i.e. [Nlon,Nlat] = size(lon)) 
% 2) om (omega) and tymes can be dimensional or nondimensional as they only appear in this functions as the product om*t.
%%

[Mm, Nn] = size(lon);
fieldMov   = zeros(Mm,Nn,length(tiltVec)); % preallocate
fieldMov_E = zeros(Mm,Nn,length(tiltVec)); % preallocate
fieldMov_N = zeros(Mm,Nn,length(tiltVec)); % preallocate

ii = 0; 

for tilt = tiltVec; ii = ii+1; tilt; display(['fraction done = ',num2str(ii/length(tiltVec))])
        fieldMov(:,:,ii) = SH2Field(Cns, lon,lat,tilt,  s,tilom,phase0, swCC); 
    if nargout > 1, 
        fieldMov_gradE(:,:,ii) =  SH2Field_gradE(Cns,    lon,lat,tilt, s,tilom,phase0, swCC); % eastward gradient
        fieldMov_gradN(:,:,ii) =  SH2Field_gradN(Cns,    lon,lat,tilt, s,tilom,phase0, swCC); % northward gradient
    end
        
end

fieldMov_gradE(:,[1 Nn],:) = 0; % correct for singularity at poles;

