clear all
close all
clc

% Given values
L = 1;                  % length of the domain
c = 1;                  % linear convection velocity
dt = 0.01;              % time step
t = 0.4;                % end time
n_t = t/dt;             % no of time steps
uspike_start = 0.1;     % value of x at which spike in velocity starts
uspike_end = 0.3;       % value of x at which spike in velocity ends
uspike_start_index = 0; % index position of node at velocity spike start
uspike_end_index = 0;   % index position of node at velocity spike end

% no of grid points
n = [20, 40, 80, 160];

% for loop to input various values of grid points
for k = 1:length(n)
    x = linspace(0, 1, n(k)); % no. of nodes
    dx = L/(n(k)-1); % node step size
    
    % for loop to find the position of node 'x' at which velocity spike
    % starts and ends
    for i = 1:n(k)
        if abs(x(i)- uspike_start) < dx/2     % using decision making statement
            uspike_start_index = i;           % -to get the position of node
        elseif abs(x(i)-uspike_end) < dx/2    % -which has the closest value to
            uspike_end_index = i;             % -to the value at x at which the
        end                                   % -velocity spike starts and ends.
    end
    u = ones(1,n(k));
    u(1) = 1;                                     % boundary condition
    u(uspike_start_index:uspike_end_index) = 2;   % assigning velocity step values
    u_old = u;
    u_initial = u;                                % initial velocity profile
    
    % for loop for time marching
    for j = 1: n_t
        % for loop for space marching
        for i = 2:n(k)
            % 1D linear convection equation, du/dt = -C(du/dx)
            u(i) = u_old(i) - (c*dt/dx)*(u_old(i)-u_old(i-1));
        end
        figure(1);
        % need fullscreen window to view 4 plots at once
        set(gcf, 'Position', get(0, 'screensize'));
        subplot(2,2,k);
        plot(x, u_initial,'r',x,u,'b');
        axis([0 1 0.5 2.5]);
        title(['No of grid points = ',num2str(n(k))]);
        xlabel('Nodes, x');
        ylabel('Velocity, u');
        legend('Original velocity profile','Final velocity profile');
        % updating u_old
        u_old = u;
        pause(0.001);
    end
end
