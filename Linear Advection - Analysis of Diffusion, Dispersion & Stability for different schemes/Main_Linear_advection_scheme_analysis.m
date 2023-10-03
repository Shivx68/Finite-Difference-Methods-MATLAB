clear all
close all
clc

% Inputs
L = 3;                  % length of the domain
v = 1;                  % linear advection velocity
t = 2;                  % end time
n = 151;                % no of grid points

x = linspace(0, L, n);  % no. of nodes
dx = L/(n-1);           % grid size
%Co = [0.1 0.2 0.4 0.8];
Co = 0.4;

for i = 1:length(Co)
  
  C = Co(i);
  dt = C*dx/v;            % time step
  n_t = t/dt;             % no of time steps

  u_initial = initialize(x);
  u_analytic = analytic(x, t, v);

  %u_FTBS(i,:) = FTBS(u_initial, n, n_t, C);
  %u_FTCS(i,:) = FTCS(u_initial, n, n_t, C);
  %u_BTCS(i,:) = BTCS(u_initial, n, n_t, C);
  %u_LAXW(i,:) = LaxWendroff(u_initial, n, n_t, C);
  %u_SOUP = SecondOrderUpwind(u_initial, n, n_t, C);
  %u_LEAP(i,:) = LeapFrog(u_initial, n, n_t, C);
  %u_RK3(i,:) = RK3UP3(u_initial, x, n_t, dt, v);
  u_RK4(i,:) = RK4UP4(u_initial, x, n_t, dt, v);
  u_RK3CD2 = RK3T_CD2S(u_initial, n, n_t, C);
  u_RK3UP2 = RK3T_UP2S(u_initial, n, n_t, C);
  %u_unstable(i,:) = RK3T_UP3S(u_initial, x, n_t, dt, v)
  u_RK4OUP4(i,:) = RK4OUP4(u_initial, x, n_t, dt, v);
  u_RK4CD4(i,:) = RK4CD4(u_initial, x, n_t, dt, v);

end

% Plotting Results
figure(1);
plot(x, u_analytic, 'k', 'LineWidth', 2);
hold on;
plot(x,u_RK4OUP4,'LineWidth', 2)%, x, u_RK4CD4,'LineWidth', 2);
%plot(x, u_RK4,'LineWidth', 2)%, x, u_SOUP,'LineWidth', 2);
for i = 1:length(Co)
  figure(1)
  plot(x, u_RK4CD4(i,:), 'LineWidth', 2);
end
hold off;
%axis([1.6 2.9 -0.5 3]);  % zoom to the pulse location at time t
%axis([0 3 -0.5 2.5]);
title('Numerical solution using RK4-OUP4 scheme');
xlabel('Domain Length');
ylabel('Scalar quantity, u');
legend('Analytic Solution', 'C = 0.1')%, 'C = 0.2', 'C = 0.4', 'C = 0.8');
%legend('Analytic Solution', 'RK4-OUP4', 'RK4-CD4', 'RK4-BUP4');
%legend('Analytic Solution', 'FTBS(1st order)', 'RK3-BUP3 (3rd order)')
%legend('Analytic Solution', 'LAXW (2nd order)', 'RK4-BUP4 (4rd order)');