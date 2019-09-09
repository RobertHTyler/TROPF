function field = SH2Field_gradN(Cns, lon,lat,t, s,om,phase0, swCC)
% function field = SH2Field_gradN(Cns, lon,lat,t, s,om,phase0, swCC)
%
% Pth (i.e. (-1)*partial_{\theta}) of mapped Cns field:
%
% IN: switchCC:      swCC = 1, unless the field expanded has a 1i*
% multiplier, in which case swCC = -1 (because we need to subtract
% rather than add the complex conjugate in the latter case)
% 
%
%% 
[Nlon,Nlat] = size(lon);
colat       = 90 - lat;     % colatitude (deg)
N           = length(Cns);  % number of SH terms in expansion
omt         = om*t;         % temporal phase omega*t



% initialize:
field      = zeros(Nlon,Nlat);
fieldslice = zeros(1,Nlat);
 
for ii=1:N;      % Sum over degrees: 
    ni=s+(ii-1); % ni is degree, because matrix sum goes ni = s, s+1,...(N-s+1)
    if ni>0 
fieldslice = fieldslice + Cns(ii)*PthPnsx(ni,s,cos(colat(1,:)*pi/180));
    end
end
 
field = repmat(fieldslice,length(Nlon),1);

phs   = s*(lon*pi/180) - omt + phase0;
field = field.*exp(1i*(phs));
field = field + swCC*conj(field);
field = -1*field; % change sign so gradient will be "North" rather than toward positive colat