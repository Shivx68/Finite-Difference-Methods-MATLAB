% Newton Raphson Function for computing M_act from f_act
function out = F(M, f, gamma)
%function out = F(M)
  gamma = 1.4;
  out = f - sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(M^2 -1))) + atand(sqrt(M^2 -1));
end
