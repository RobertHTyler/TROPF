function ratiofactorials1 = ratiofactorials1(n,s);
% function ratiofactorials1 = ratiofactorials1(n,s);
%
% Function to stably calculate ratio (n+s)!/(n-s)! 
% 
% In: n (degree), s (order/rank)
% Out: ratio (n+s)!/(n-s)! 
%
% Notes:
% 1) Matlab calls the factorial operator ! by name "factorial".
% 
% R. Tyler (10 Aug. 2017)
% %

% initialize:
w = 0; produkt = 1;

for ii = 0:2*s-1
    junk = s - ii;
    fac  = (n+junk);
    produkt = produkt.*fac;
end

ratiofactorials1 = produkt;



