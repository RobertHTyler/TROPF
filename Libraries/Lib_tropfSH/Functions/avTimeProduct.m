function   field_avTimeProduct = avTimeProduct(cns1,s1,om1,   cns2,s2,om2, lon,lat) 
% function   field_avTimeProduct = avTimeProduct(cns1,s1,om1,   cns2,s2,om2, lon,lat) 
%
% 
%%
% 
% cns1 = -1i*Gns; 
% s1   = sF;
% om1  = tilom;
% 
% 
% cns2 = -1i*Gns;
% s2   = sF;
% om2  = tilom;
% 


[Nlon,Nlat] = size(lon); % grid dimensions

N = length(cns1); % number of SH terms in expansion
omjunk     = 0; % this is multiplied by t which is zero anyway
phase0junk = 0; 

% field 1:
var_Ynsom_arr = Ynsom_arr(lon,lat,0, N, s1,omjunk,phase0junk);
field1        = permute( reshape( var_Ynsom_arr*cns1 , Nlat,Nlon ), [2 1] ); % complex field from SH coefs

% field 2:
var_Ynsom_arr = Ynsom_arr(lon,lat,0, N, s2,omjunk,phase0junk);
field2        = permute( reshape( var_Ynsom_arr*cns2 , Nlat,Nlon ), [2 1] ); % complex field from SH coefs

if om1 == om2,
fac1 = 0; fac2 = 1;
elseif om1 == -om2,
fac1 = 1; fac2 = 0;
else 
fac1 = (1i/(2*pi))*(exp(-1i*2*pi*om2/om1)-1)/(1+om2/om1+eps);
fac2 = (1i/(2*pi))*(exp( 1i*2*pi*om2/om1)-1)/(1-om2/om1+eps);
end


junk  = (field1.*field2) * fac1   + (field1.*conj(field2)) * fac2 ;
  
 
field_avTimeProduct = junk + conj(junk);
 

