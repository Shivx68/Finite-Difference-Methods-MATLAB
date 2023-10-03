function [u] = FTBS(u_initial, n, n_t, C)
% Forward Time Backward Space - 1st order in time and space
  
  u = u_initial;
  % for loop for time marching
  for j = 1: n_t
    % updating u_old
    u_old = u;
    % for loop for space marching
    for i = 2:n
        % 1D linear convection equation, du/dt = -C(du/dx)
        u(i) = u_old(i) - C*(u_old(i)-u_old(i-1));
    end
    
    % Boundary condition
    u(1) = 1;    
    
  end
end