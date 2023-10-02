function [ydot] = odefun_pendulum(t,y,L,m,b,g)
% defining size of ydot matrix
ydot = zeros(2,1);
% d(theta)/dt = omega (angular velocity)
ydot(1) = y(2);
% d(omega) = alpha = -b*omega/m - g*sin(theta)/L
ydot(2) = -(b*y(2))/m - (g*sin(y(1)))/L;
end