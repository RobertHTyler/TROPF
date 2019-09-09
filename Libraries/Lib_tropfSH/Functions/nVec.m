function [nvec,Ntrunc] = nVec(N,s);
% function [nvec,Ntrunc] = nVec(N,s);
% Returns a column vector of the SH degrees involved in an expansion with N terms
%
% Input: 
% N       : Number of SH terms in the expansion (scalar)
% s       : Order/rank of expansion terms (scalar)
%
% Output:
% nvec    : Column vector of the SH degrees spanning  n = s, s+1, s+2...Ntrunc
% Ntrunc  : Truncation degree
% 
% R. Tyler, 28 Dec., 2018

Ntrunc  = N + s - 1         ; % truncation degree of SH expansion
nvec    = [s:Ntrunc]'       ; % vector of SH degrees involved
