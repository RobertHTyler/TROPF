function mabs = mabs(A1,A2)
% Returns the max of the absolute value of differences between arrays A1,
% A2
%
% 
mabs = max(abs(A1(:)-A2(:)));