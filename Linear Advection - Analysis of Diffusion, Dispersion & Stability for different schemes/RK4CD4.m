function [u] = RK4CD4(u_initial, x, n_t, dt, v)
  
  u = u_initial;
  n = length(x);
  dx = x(2)-x(1);
  
  % for loop for time marching
  for j = 1: n_t
    
    % updating u_old
    u_old = u;
    
    % Runge Kutta- First Step
    for i = 3:n-2
        % 1D linear convection equation, du/dt = -C(du/dx)
        k1(i) = (-v/(12*dx))*(-u_old(i+2)+8*u_old(i+1)-8*u_old(i-1)+u_old(i-2));
        u1(i) = u_old(i) + (k1(i)/2)*dt;
    end
    % Boundary condition
    u1(1) = 1;
    u1(2) = 1;
    u1(n-1) = 1;  
    u1(n) = 1;
    
    % Runge Kutta- Second Step
    for i = 3:n-2
        % 1D linear convection equation, du/dt = -C(du/dx)
        k2(i) = (-v/(12*dx))*(-u1(i+2)+8*u1(i+1)-8*u1(i-1)+u1(i-2));
        u2(i) = u_old(i) + (k2(i)/2)*dt;
    end
    % Boundary condition
    u2(1) = 1;
    u2(2) = 1;
    u2(n-1) = 1;  
    u2(n) = 1;
    
    % Runge Kutta- Third Step
    for i = 3:n-2
        % 1D linear convection equation, du/dt = -C(du/dx)
        k3(i) = (-v/(12*dx))*(-u2(i+2)+8*u2(i+1)-8*u2(i-1)+u2(i-2));
        u3(i) = u_old(i) + k3(i)*dt;
    end
    % Boundary condition
    u3(1) = 1;
    u3(2) = 1;
    u3(n-1) = 1;  
    u3(n) = 1;
    
    % Runge Kutta- Fourth Step
    for i = 3:n-2
        % 1D linear convection equation, du/dt = -C(du/dx)
        k4(i) = (-v/(12*dx))*(-u3(i+2)+8*u3(i+1)-8*u3(i-1)+u3(i-2));
        u(i) = u_old(i) + (dt/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
    end
    % Boundary condition
    u(1) = 1;
    u(2) = 1;
    u(n-1) = 1;   
    u(n) = 1;

  end
end