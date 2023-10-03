function [u] = LeapFrog(u_initial, n, n_t, C)
% LeapFrog scheme- Second Order Central difference in both time and space
  
  u = u_initial;
  u_old = u;
  
  % First step initialization using upwind
  % for loop for space marching
  for i = 2:n
      % 1D linear convection equation, du/dt = -C(du/dx)
      u(i) = u_old(i) - C*(u_old(i)-u_old(i-1));
  end
  u(1) = u(2);
  % Initializing for 2 time levels
  u_older = u_old;
  u_old = u;
  
  % for loop for time marching
  for j = 1: n_t
    % updating u_old and u_older
    u_older = u_old;
    u_old = u;
    
    % for loop for space marching
    for i = 2:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        u(i) = u_older(i) - C*(u_old(i+1)-u_old(i-1));
    end
    % Boundary condition
    u(1) = 1;
    u(n) = 1;
    
  end
end
