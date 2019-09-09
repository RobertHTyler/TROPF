function [Dxvar,Dyvar,Gradstruct] = gradGlobe(lon,lat,radius,var)
%% 
% Dxvar,Dyvar : gradient at center points (step over two grid spaces)
% Gradstruct  : gradient at faces (step over one grid space)
% Assumes regular grid (global lon,lat (degrees)) with no overlap points in E/W boundary
%
% Determine resolution:
dellon = lon(2,1)-lon(1,1); 
dellat = lat(1,2)-lat(1,1);

% Make overlapping vars/params:

lonext = [lon(1,:)-dellon; lon; lon(end,:)+dellon]; % lon extended in lon direction
latext4lon = [lat(1,:); lat; lat(end,:)];    % lat extended in lon direction
latext = [lat(:,1)-dellat, lat, lat(:,end)+dellat]; % lat extending in lat direction
%
varlonext = [var(end,:); var; var(1,:)];
varlatext = [var(:,1), var, var(:,end)];
%
x=radius*(lonext*pi/180).*cos(latext4lon*pi/180);
y=radius*(latext*pi/180);

% Calc gradients:
Dxvar = (varlonext(3:end,:) - varlonext(1:end-2,:)) ./ (x(3:end,:) - x(1:end-2,:)); 
Dyvar = (varlatext(:,3:end) - varlatext(:,1:end-2)) ./ (y(:,3:end) - y(:,1:end-2));


% Calc gradients:
Gradstruct.w = (varlonext(2:end-1,:) - varlonext(1:end-2,:)) ./ (x(2:end-1,:) - x(1:end-2,:)); 
Gradstruct.e = (varlonext(3:end,:) - varlonext(2:end-1,:)) ./ (x(3:end,:) - x(2:end-1,:)); 
Gradstruct.s = (varlatext(:,2:end-1) - varlatext(:,1:end-2)) ./ (y(:,2:end-1) - y(:,1:end-2));
Gradstruct.n = (varlatext(:,3:end) - varlatext(:,2:end-1)) ./ (y(:,3:end) - y(:,2:end-1));
