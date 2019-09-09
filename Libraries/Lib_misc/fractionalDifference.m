function fractionalDifference = fractionalDifference(v1,v2);
% Returns the fractional difference between two vectors
fractionalDifference = max(abs(v1-v2)) / mean( 0.5*(abs(v1) + abs(v2)) );