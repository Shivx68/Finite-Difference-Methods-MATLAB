function [u] = RK4OUP4(u_initial, x, n_t, dt, v)
  
  u = u_initial;
  n = length(x);
  dx = x(2)-x(1);
  
  % for loop for time marching
  for j = 1: n_t
    
    % updating u_old
    u_old = u;
    
    % Runge Kutta- First Step
    for i = 5:n
        % 1D linear convection equation, du/dt = -C(du/dx)
        k1(i) = (-v/(12*dx))*(25*u_old(i)-48*u_old(i-1)+36*u_old(i-2)-16*u_old(i-3)+3*u_old(i-4));
        u1(i) = u_old(i) + (k1(i)/2)*dt;
    end
    % Boundary condition
    u1(1) = 1;
    u1(2) = 1;
    u1(3) = 1;  
    u1(5) = 1;
    
    % Runge Kutta- Second Step
    for i = 5:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k2(i) = (-v/(12*dx))*(25*u1(i)-48*u1(i-1)+36*u1(i-2)-16*u1(i-3)+3*u1(i-4));
        u2(i) = u_old(i) + (k2(i)/2)*dt;
    end
    % Boundary condition
    u2(1) = 1;
    u2(2) = 1;
    u2(3) = 1;  
    u2(4) = 1;
    
    % Runge Kutta- Third Step
    for i = 5:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k3(i) = (-v/(12*dx))*(25*u2(i)-48*u2(i-1)+36*u2(i-2)-16*u2(i-3)+3*u2(i-4));
        u3(i) = u_old(i) + k3(i)*dt;
    end
    % Boundary condition
    u3(1) = 1;
    u3(2) = 1;
    u3(3) = 1;  
    u3(4) = 1;
    
    % Runge Kutta- Fourth Step
    for i = 5:n-1
        % 1D linear convection equation, du/dt = -C(du/dx)
        k4(i) = (-v/(12*dx))*(25*u3(i)-48*u3(i-1)+36*u3(i-2)-16*u3(i-3)+3*u3(i-4));
        u(i) = u_old(i) + (dt/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
    end
    % Boundary condition
    u(1) = 1;
    u(2) = 1;
    u(3) = 1;   
    u(4) = 1;

  end
end