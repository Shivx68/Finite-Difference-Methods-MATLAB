% Function to solve 2D Transient state heat conduction equation explicitly

function transient_explicit(L, n, t, dt)

    x = linspace(0,L,n);        % x nodes
    y = linspace(0,L,n);        % y nodes
    dx = L/(n-1);               % gird size along x
    dy = L/(n-1);               % grid size along y
    n_t = t/dt;                 % number of time steps
    %alpha = 1;                  % thermal diffusivity 
    alpha = 0.045

    %Initialization
    T = 300*ones(n, n);             % initializing T matrix
    T(:,1) = 400;               % left boundary condition
    T(:,n) = 800;             % right boundary condition
    T(1,:) = 900;               % bottom boundary condition
    T(n,:) = 600;             % top boundary condition
    T_old = T;                  % for updation old values in convergence loop

    k1 = alpha*dt/dx^2;         % for ease of calculation
    k2 = alpha*dt/dy^2;

    % time loop    
    for k = 1:5000 %n_t   
        % nodal loop
        for j = 2: (n-1)
            for i = 2: (n-1)
                % explicit scheme for 2D heat equation
                T(i,j) = (1-2*k1-2*k2)*T_old(i,j)+k1*(T_old(i+1,j)+T_old(i-1,j))...
                    +k2*(T_old(i,j+1)+T_old(i,j-1));
            end
        end
        % updating old values
        T_old = T;
        
        %pause(0.0001)
    end
    % creating a contour plot
        [C,h] = contourf(x,y,T);
        colormap(jet);
        clabel(C,h);
        xlabel('X axis');
        ylabel('Y axis');
        title_text = sprintf(['Solving transient 2D heat equation explicitly'...
            '\nNo. of grid points = %d; dx = %0.5g; dt = %g'...
            '\nTemp. distribution at time, t= %g'], n, dx, dt, dt*k);
        title(title_text);
end