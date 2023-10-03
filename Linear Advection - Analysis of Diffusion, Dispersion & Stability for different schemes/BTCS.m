function [u] = BTCS(u_initial, n, n_t, C)
% Backward(Implicit) Time Central Space - 1st order time, 2nd order space

  u = u_initial;          % initial velocity profile
  tol = 1e-4;             % tolerance for convergence
  % for loop for time marching
  for j = 1: n_t
    % updating prev time step u
    u_old = u;
    error = 1;
    % Convergence loop
    while(error>tol)
      % updating last convergence iteration u
      u_guess = u;
      % for loop for space marching
      for i = 2:n-1
          % 1D linear convection equation, du/dt = -C(du/dx)
          u(i) = u_old(i) - (0.5*C)*(u_guess(i+1)-u_guess(i-1));  % Gauss-Seidel solver
      end
      % computing error
      error = max(abs(u-u_guess))
    end
    disp(j);
    % Boundary condition
    u(1) = 1;
    u(n) = 1;
   
  end
end