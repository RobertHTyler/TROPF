
% This script sets the TROPF paths and includes any general startup and
% housekeeping tasks. 
%
% Instructions for beginners:
%
% Define TROPFPATH below (usually you can set this to be your path to the
% directory in which this startup script is located). The other directories
% will be automatically set based on this.
%
%
% R. Tyler (1 Nov. 2018)


% % Define path to TROPF root directory:
TROPFPATH = '~/Desktop/TROPF'; 


% % Add other paths to TROPF libraries:
%
% add path to TROPF Spherical-Harmonic library:
addpath([TROPFPATH,'/Libraries/Lib_tropfSH/Functions']); 
addpath([TROPFPATH,'/Libraries/Lib_tropfSH/Scripts']); 
%
% add path to TROPF Finite-Volume library:
addpath([TROPFPATH,'/Libraries/Lib_tropfFV']);
%
% add path to miscellaneous library:
addpath([TROPFPATH,'/Libraries/Lib_misc']);   
addpath([TROPFPATH,'/Libraries/Lib_misc/othercolor']);  % colormaps from U. Oregon (see notes below and license.txt)
addpath([TROPFPATH,'/Libraries/Lib_misc/Colormaps']);   % colormaps from Python forums (see notes below and license.txt)




% Notes on othercolor library:
%
% 1) Run colormap(othercolor('BuDRd_18')) for nice alternative to 'jet'
% when you want to highlight extrema in data 
%
% 2) This is obtained from
% http://geography.uoregon.edu/datagraphics/color_scales.htm  


% Notes on Colormaps library:
%
% 1) Obtained from https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps
