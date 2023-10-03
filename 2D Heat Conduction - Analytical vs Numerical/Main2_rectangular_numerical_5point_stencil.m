% 3A - Solving the Laplacian (Temperature steady-state) using 5-point stencil
clear all;
close all;
clc;

% Matrix Dimensions
L = 100;
H = 100;
T = zeros(L,H);

% Boundaries
T(1,:) = 0;
T(L, :) = 0;
T(:,1) = 0;
T(:, H) = 100;

n = 1;
% for loop to reach convergence 
for k = 1:1000
  
  % Nodal loop
  for i = 2:H-1
    for j = 2:L-1
      T(i,j) = (1/4)*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
    end
  end
  
  % Contour Plot
  contourf(T,100, 'LineStyle', 'none');
  colormap(jet);
  xlabel('X axis');
  ylabel('Y axis');
  title_text = sprintf('Numerical solution using 5-point stencil, n-iter = %d', n);
  title(title_text);
  drawnow();
  n = n+1;
end

