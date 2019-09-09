function LVi = build_LVi_fromSlowness(N, tilnusqns);
%% function LVi = build_LVi_fromSlowness(N, tilnusqns);
% Builds the matrix LVi, which represents the operator ($ L_Vi $), from
% slowness parameter vector tilnsqns ($ (\tilde{\nu}_s^n)^2 $ )
% (tilnusqns may be scalar or column vector; tilnusqns may be complex)
%
% Input: 
% tilnusqns : Column vector of the slowness parameter (scalar or column vector; may be complex)
%             If tilnusqns is a scalar, LVi is a diagonal matrix repeating this scalar value
%             If tilnusqns is a vector, LVi is a diagonal matrix using this vector (i.e. slowness varies with degree)
%        
%            
%
% Output:
% LVi       : Square matrix (sparse, diagonal)
% 
% R. Tyler, 28 Dec., 2018
%%

if length(tilnusqns)>1
LVi  = spdiags( 1./tilnusqns  , 0,N,N) ;  
else
LVi  = spdiags( 1./(tilnusqns*ones(N,1)) , 0,N,N) ;
end
