function varNotStruct = loadVar(fylename,var)
% load variable as array rather than structure
% (when Matlab load is used with an output argument, it loads 
% (even if only one array is requested) into a structure. This program just
% converts it back into an array in the workspace)
%
% R. Tyler 1 Aug. 2018


structNotVar = load(fylename,var);
varNotStruct = struct2cell(structNotVar);
varNotStruct = varNotStruct{1};
