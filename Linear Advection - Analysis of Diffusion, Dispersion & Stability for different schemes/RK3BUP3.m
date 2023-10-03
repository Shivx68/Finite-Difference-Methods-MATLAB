   function [u] = RK3BUP3(u_initial, x, n_t, dt, v)

  u = u_initial;
  n = length(x);
  dx = x(2)-x(1);
  
  % for loop for time marching
  for j = 1: n_t
    % updating u_old
    u_old = u;
    
    % Runge Kutta- First Step
    for i = 3:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k1(i) = (-v/(6*dx))*(2*u_old(i+1)+3*u_old(i)-6*u_old(i-1)+u_old(i-2));
        u1(i) = u_old(i) + (k1(i)/2)*dt;
    end
    % Boundary condition
    u1(1) = 1;
    u1(2) = 1; 
    u1(n) = 1;
    % Runge Kutta- Second Step
    for i = 3:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k2(i) = (-v/(6*dx))*(2*u1(i+1)+3*u1(i)-6*u1(i-1)+u1(i-2));
        u2(i) = u_old(i) - k1(i)*dt + 2*k2(i)*dt;
    end
    % Boundary condition
    u2(1) = 1;
    u2(2) = 1;  
    u2(n) = 1;
    % Runge Kutta- Third Step
    for i = 3:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k3(i) = (-v/(6*dx))*(2*u2(i+1)+3*u2(i)-6*u2(i-1)+u2(i-2));
        u(i) = u_old(i) + (dt/6)*(k1(i)+4*k2(i)+k3(i));
    end
    % Boundary condition
    u(1) = 1;
    u(2) = 1;   
    u(n) = 1;
  end
end