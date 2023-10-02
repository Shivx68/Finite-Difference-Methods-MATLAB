% Function to solve 2D Transient state heat conduction implicitly using Jacobi
% iterative solver

function transient_implicit_jacobi(L, n, tolerance, t, dt)

    x = linspace(0,L,n);        % x nodes
    y = linspace(0,L,n);        % y nodes
    dx = L/(n-1);               % gird size along x
    dy = L/(n-1);               % grid size along y
    n_t = t/dt;                 % number of time steps
    alpha = 1;                  % thermal diffusivity 


    %Initialization
    T = ones(n, n);             % initializing T matrix
    T(:,1) = 400;               % left boundary condition
    T(:,n) = 800;               % right boundary condition
    T(1,:) = 900;               % bottom boundary condition
    T(n,:) = 600;               % top boundary condition
    T_old = T;                  % for updation old values in convergence loop
    T_prev = T;
    n_iteration = 1;            % to count the number of total iterations            

    k1 = alpha*dt/dx^2;         % for ease of calculation
    k2 = alpha*dt/dy^2;
    
    % time loop    
    for k = 1:n_t
        error = 1;              % error initialized to 1 before each time loop
        % convergence loop
        while error > tolerance
            % nodal loop
            for j = 2: (n-1)
                for i = 2: (n-1)
                    % implicit scheme using Jacobi
                    T(i,j) = (T_prev(i,j)+k1*(T_old(i-1,j)+T_old(i+1,j))+...
                        k2*(T_old(i,j-1)+T_old(i,j+1)))/(1+2*k1+2*k2);
                end
            end
            % convergence criterion
            error = max(max(abs(T - T_old)));
            % updating old values
            T_old = T;
            n_iteration = n_iteration +1;
        end
        % updating previous values
        T_prev = T;
        % creating a contour plot
        [C,h] = contourf(x,y,T);
        clabel(C,h);
        colormap(jet);
        xlabel('X axis');
        ylabel('Y axis');
        title_text = sprintf(['Solving transient 2D heat equation implicitly'...
            ' using Jacobi\nNo. of grid points = %d; dx = %0.5g; dt = %g'...
            '\nTemp. distribution at time, t= %g; Total iterations= %g']...
            , n, dx, dt, dt*k, n_iteration);
        title(title_text);
        pause(0.0001)
    end
end