function [tilcesq,tilalpp] = cNalppFromSlow(tilom,tilnusqns)
%
% Based on definition in tilnusqns =  (1 + 1i*tilalpp/tilom)./tilcesq;
% build_LV_fromSlowness


tilcesq =  1./real(tilnusqns) ;
tilalpp =  tilom * imag(tilnusqns)/real(tilnusqns);



