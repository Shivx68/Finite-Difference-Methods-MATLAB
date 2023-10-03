% 2A - Analytical Solution using a Convenient Boundary Condition
clear all;
close all;
clc;

% Matrix dimensions (Domain)
L = 100;
H = 100;

% Declaring T matrix
T = zeros(L,H);

% Grid spacing along x and y
dx = 1;
dy = 1;

% Defining the no of points along x and y
x = 1:dx:L;
y = 1:dy:H;

% Obtain the X and Y locations for each node
[X,Y] = meshgrid(x, y);

% Boundary Conditions
BC = sin(2*pi*y/H);

% Solution
A2 = 2/(H*sinh(2*pi*L/H))*sum(BC.^2);
T = A2*sinh(2*pi*X/H).*sin(2*pi*Y/H);

% Contour Plot
[C,h] = contourf(T,100, 'LineStyle', 'none');
colormap(jet);
xlabel('X axis');
ylabel('Y axis');
title_text = sprintf('Analytical(Infinite-time) Solution- easy boundary condition');
title(title_text);
