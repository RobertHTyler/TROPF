function [lowval,highval,rangevals] = rangelh(var);
% extrema and range in whole data set var:
lowval    = min(var(~isnan(var))), 
highval   = max(var(~isnan(var))), 
rangevals = highval-lowval,