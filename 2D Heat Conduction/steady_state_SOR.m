% Function to solve 2D Steady state heat conduction equation using Successive Over
% Relaxation

function steady_state_SOR(L, n, tolerance, omega)
    
    % Constants
    x = linspace(0,L,n);            % x nodes
    y = linspace(0,L,n);            % y nodes
    dx = L/(n-1);                   % gird size along x
    dy = L/(n-1);                   % grid size along y

    %Initialization
    T = ones(n, n);                 % initializing T matrix
    T(:,1) = 400;                   % left boundary condition
    T(:,n) = 800;                   % right boundary condition
    T(1,:) = 900;                   % bottom boundary condition
    T(n,:) = 600;                   % top boundary condition
    T_old = T;                      % for updation old values in convergence loop
    error = 1;                      % to enter convergence loop first time
    n_iteration = 1;                % to record no. of iterations
    k = 2*(dx^2+dy^2)/(dx^2*dy^2);  % for ease of calculation

    T_gs = zeros(n, n);         

     %convergence loop
    while error > tolerance
        % nodal loop
        for j = 2: (n-1)
            for i = 2: (n-1)
                % sor method- advanced gauss-seidel method
                T_gs(i,j) = (1/k)*((T(i-1,j)+T_old(i+1,j))/dx^2+...
                    (T(i,j-1)+T_old(i,j+1))/dy^2);
                T(i,j) = T_old(i,j)*(1-omega)+omega*T_gs(i,j);
            end
        end
        % checking convergence
        error = max(max(abs(T - T_old)));
        %updating old values
        T_old = T;
        n_iteration = n_iteration + 1;
        
        
    end
    % creating a contour plot
        [C,h]=contourf(x,y,T,100,'LineStyle', 'none');
        colormap(jet);
        %colorbar;
        %clabel(C,h);
        %xlabel('X axis');
        %ylabel('Y axis');
        %title_text = sprintf(['Solving 2D steady state heat equation using'...
            %' SOR iterative solver\nNo. of grid points = %d;'...
            %' Total iterations = %d'], n, n_iteration);
        %title('Temperature at Steady-State, K');
        pause(0.0001);
    
end
% iterations number increases above and below omega = 1.73 this value and after 1.9 the solution does not converge
% for higher no. of grid points the optimum omega value increases & no. of
% iterations also increases
% when omega >= 2, the solution does not converge