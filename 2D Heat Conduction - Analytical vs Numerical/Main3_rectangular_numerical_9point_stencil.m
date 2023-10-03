% 3B- Solving the Laplacian (Temperature steady-state) using 9-point stencil
clear all;
close all;
clc;

% Creating a 2D rectangular domain

L = 100; % no of nodes along y
H = 100; % no of nodes along x
T = zeros(L,H);

% Applying the Boundary Conditions

T(1,:) = 0;   % top face
T(L, :) = 0;  % bottom face
T(:,1) = 0;   % left face
T(:, H) = 1;  % right face

n_iter = 1000; % no. of iterations
n = 1;         % iteration counter
% for loop to reach convergence 
for n = 1:n_iter
  % Computing the solution for each node
  for i = 2:H-1
    for j = 2:L-1
      T(i,j) = (1/8)*(T(i+1,j) + T(i-1,j) + T(i,j+1)...
      + T(i,j-1) + T(i+1,j+1) + T(i-1,j+1) + T(i+1, j-1) + T(i-1, j-1));
    end
  end
  % Contour plot
  [C,h] = contourf(T,100, 'Linestyle', 'none');
  colormap(jet);
  xlabel('X axis');
  ylabel('Y axis');
  title_text = sprintf('Numerical solution using 9-point stencil, n-iter = %d', n);
  title(title_text);
  drawnow();
  n = n+1;
end
