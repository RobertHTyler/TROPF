function field = SH2field_gradE(Cns, lon,lat,t, s,om,phase0, swCC)
% function field = SH2field_gradE(Cns, lon,lat,t, s,om,phase0, swCC)
% Returns Pph (i.e. (1/sin(\theta)*partial_{\phi}) of mapped Cns field:
%
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
fieldslice = fieldslice + 1i*Cns(ii)*sPnsxOVsinth(ni,s,cos(colat(1,:)*pi/180));
    end
end

field = repmat(fieldslice,length(Nlon),1);

phs   = s*(lon*pi/180) - omt + phase0;
field = field.*exp(1i*(phs));
field = field + swCC*conj(field);
field(:,[1,end]) = real(field(:,[1,end])); % get rid of imaginary part of NaN at poles (real part of NaN remains)

