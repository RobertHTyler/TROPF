function LV = build_LV_fromSlowness(N,tilnusqns);
%% function LV = build_LV_fromSlowness(N,tilnusqns);
% Builds the matrix LV, which represents the operator ($ L_V $), from
% slowness parameter vector tilnsqns ($ (\tilde{\nu}_s^n)^2 $ )
% (tilnusqns may be scalar or column vector; tilnusqns may be complex)
%
% Input: 
% tilnusqns : Column vector of the slowness parameter (scalar or column vector; may be complex)
%             If tilnusqns is a scalar, LV is a diagonal matrix repeating this scalar value
%             If tilnusqns is a vector, LV is a diagonal matrix using this vector (i.e. slowness varies with degree)
%        
%            
%
% Output:
% LBi     : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018
%%

if length(tilnusqns)>1
LV  = spdiags( tilnusqns  , 0,N,N) ;  
else
LV  = spdiags( tilnusqns*ones(N,1) , 0,N,N) ;
end
