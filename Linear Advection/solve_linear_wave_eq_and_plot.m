% Function that takes time step sizes as an argument to slove the linear
% convection equation and plots the velocity profile for each value of step
% size

function solve_linear_wave_eq_and_plot(dt)

    % Given values
    L = 1;                  % length of the domain
    c = 1;                  % linear convection velocity
    spike_start = 0.1;      % value of x at which spike in velocity starts
    spike_end = 0.3;        % value of x at which spike in velocity ends
    n = 80;                 % no. of grid points

    x = linspace(0, 1, n);  % no. of nodes
    dx = L/(n-1);           % node step size

    % calling function to find the index positions where velocity step starts and ends
    spike_start_index = find_index_position(x, spike_start);
    spike_end_index = find_index_position(x, spike_end);

    t = 0.4;                         % end time

    % for loop to input various values of time step size
    for k = 1:length(dt)
    
        n_t = t/dt(k);               % no of time steps
        u = ones(1, n);
        u(1) = 1;
        u(spike_start_index:spike_end_index) = 2;
        u_initial = u;
        u_old = u;
        % for loop for time marching
        for j = 1: n_t

            % for loop for space marching
            for i = 2:n
                % 1D linear convection equation, du/dt = -c(du/dx)
                u(i) = u_old(i) - (c*dt(k)/dx)*(u_old(i)-u_old(i-1));
            end
            
            % updating u_old
            u_old = u;
        end
        figure(1);
        % need fullscreen window to view 4 plots at once
        set(gcf, 'Position', get(0, 'screensize'));
        subplot(2,2,k);
        plot(x, u_initial,'r',x,u,'b');
        axis([0 1 0.5 2.5]);
        title(['Courant No = ',num2str(c*dt(k)/dx)]);
        xlabel('Nodes, x');
        ylabel('Velocity, u');
        legend('Original velocity profile','Final velocity profile');
    end
