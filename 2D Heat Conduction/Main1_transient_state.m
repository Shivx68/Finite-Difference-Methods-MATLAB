% Solving 2D transient state heat conduction equation using explicitly and 
% implicitly using Jacobi, Gauss-Seidel, and SOR iterative solvers

clear all
close all
clc

L = 1;                      % domain length along x and y                                          
%n = 20;                     % number of grid points along x and y
n = 6;
tolerance = 1e-4;           % convergence criterion
t = 2;                  % desired time of temperature distribution
%dt = 1e-4;                  % time step size
dt = 0.25;
omega = 1.8;                % relaxation factor
% method = 1, for solving Explicitly
% method = 2, for solving Implicitly using Jacobi
% method = 3, for solving Implicitly using Gauss-Seidel
% method = 4, for solving Implicitly using SOR

method = 4;
omega = linspace(1.23,1.25,20)
%omega = 1.7;
for i = 1:length(omega)
  
  if method == 1
      transient_explicit(L, n, t, dt); %function for solving using explicit scheme
  elseif method == 2
      transient_implicit_jacobi(L, n, tolerance, t, dt); % implicit jacobi
  elseif method == 3
      transient_implicit_gauss_seidel(L, n, tolerance, t, dt); % implicit gauss-seidel
  elseif method == 4
      iter(i) = transient_implicit_SOR(L, n, tolerance, t, dt, omega(i)); % implicit SOR
  elseif method == 5
      iter(i) = transient_implicit_SOR_jac(L, n, tolerance, t, dt, omega(i)); % implicit SOR-jacobi
  end
  disp(omega(i));
endfor
plot(omega, iter);