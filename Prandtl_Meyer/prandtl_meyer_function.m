% Prandtl-Meyer Function

function f = prandtl_meyer_function(M,gamma)
  f = sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1))) - atan(sqrt(M^2-1));
end
