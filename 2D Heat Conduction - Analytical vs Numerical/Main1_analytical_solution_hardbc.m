% 2B - Analytical Solution using a Hard Boundary Condition
clear all;
close all;
clc;

% Matrix Dimensions
L = 100;
H = 100;
T = zeros(L,H);

[X,Y] = meshgrid(1:1:L, 1:1:H);

% Boundary Condition
BC = 100;

% Solution
for k = 1:99
  Ak = 2/(H*sinh(k*pi*L/H))*sum(BC.*sin(k*pi*(1:H)/H));
  T = T + Ak*sinh(k*pi*X/H).*sin(k*pi*Y/H);
  
end

% Contour Plot
contourf(T,100, 'LineStyle', 'none');
colormap(jet);
xlabel('X axis');
ylabel('Y axis');
title_text = sprintf('Analytical(Infinite-time) Solution- hard boundary condition');
title(title_text);