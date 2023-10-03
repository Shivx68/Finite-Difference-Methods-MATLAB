function [u] = SecondOrderUpwind(u_initial, n, n_t, C)
%One sided upwind difference, 2nd order in both time and space

  u = u_initial;
  % for loop for time marching
  for j = 1: n_t
    % updating u_old
    u_old = u;
    % for loop for space marching
    for i = 4:n
        % 1D linear convection equation, du/dt = -C(du/dx)
        u(i) = u_old(i) - (C/2)*(3*u_old(i)-4*u_old(i-1)+u_old(i-2))...
          + (C^2/2)*(2*u_old(i)-5*u_old(i-1)+4*u_old(i-2)-u_old(i-3));
    end
    
    % Boundary condition
    u(1) = 1;    

  end
end