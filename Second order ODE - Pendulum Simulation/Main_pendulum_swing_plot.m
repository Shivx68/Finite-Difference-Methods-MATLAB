clear all
close all
clc

% Given values
L = 1;
m = 1;
b = 0.05;
g = 9.81;

t = [0, 20]; % time span (0 - 20) seconds
y0 = [0, 3]; % initial conditions for angular displacement and velocity

% calling function to solve the governing equation
[t, y] = ode45(@(t,y)odefun_pendulum(t,y,L,m,b,g),t,y0); 

% creating a file object to record the animation in a new file
mov = VideoWriter('Pendulum_animation.avi');
open(mov); % opening the file object to allow recording

figure(1)
% creating objects for fixed point, arm and bob of the pendulum so that
% they can be called later for updating their co-ordinates
% plotting the fixed point on the origin
fxd_point = plot(0,0,'-s','MarkerSize',5,'LineWidth',5,'Color','red');
hold on
% plotting the line from the origin to the Bob
arm = plot([0,0], [0,-1],'LineWidth',2,'Color','blue');
hold on
% plotting Bob marker
bob = plot(0,-1,'o','MarkerSize',10,'LineWidth',10,'Color','black');
axis([-2 2 -2 2]); % defining axes size
grid on; % enabling grid

% loop for updating the co-ordinates of the above plots to create an animation    
for i = 1:length(t)
    % the x and y coordinates of bob change with angular displacement as
    % given in the formula below
    x1= L*sin(y(i,1));
    y1= -L*cos(y(i,1));
    
    % updating the coordinates data for arm and bob
    set(arm,'XData',[0,x1], 'YData',[0,y1]);
    set(bob,'XData',x1, 'YData',y1);
    drawnow % command to update the output window with each iteration
    
    frame = getframe(gcf); % updated figure captured on each iteration
    writeVideo(mov, frame); % store each frame inside the object as a movie
end
close(mov); % closing the file object to stop recording

