function ST_globeTimeAverage = globeTimeAverage(Ans,Bns, s);
% function ST_globeTimeAverage = globeTimeAverage(Ans,Bns, s);
%
% Returns globe/time average of product of real fields 
% S and T (with associated coefs Ans, Bns) having
% common s and omega (degree and frequency).
%
% Input: 
% Ans:  SH coefs for real field S
% Bns:  SH coefs for real field T
%
% Output:
% ST_globeTimeAverage: 
% 
% R. Tyler (28 Dec. 2018)

N       = length(Ans) ; % number of terms in SH expansion
Ntrunc  = N+s-1       ; % truncation degree of SH expansion
nvec    = [s:Ntrunc]' ; % vector of SH degrees involved

ST_globeTimeAverage = ( (Ans).*conj(Bns) + conj(Ans).*(Bns) ) .*  (1./(2*nvec+1)).*ratiofactorials1(nvec,s);


