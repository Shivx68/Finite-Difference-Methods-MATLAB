function [u] = FTCS(u_initial, n, n_t, C) 
% Forward Time Central Space - 1st order in time, 2nd order in space

  u = u_initial;
  % for loop for time marching
  for j = 1: n_t
    % updating u_old
    u_old = u;
    % for loop for space marching
    for i = 2:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        u(i) = u_old(i) - (0.5*C)*(u_old(i+1)-u_old(i-1));
    end
    % Boundary condition
    u(1) = 1;    
    u(n) = 1;
  
  end
end