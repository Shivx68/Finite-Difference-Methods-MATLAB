% Solving steady state 2D heat conduction equation over a square domain
% using Jacobi, Gauss_Seidel and SOR iterative solvers

clear all
close all
clc

L = 1;                      % domain length along x and y                                          
n = 20;                     % number of grid points along x and y
tolerance = 1e-4;           % convergence criterion
omega = 1.73;                % over relaxation factor

% method = 1, for Jacobi
% method = 2, for Gauss-Seidel
% method = 3, for SOR
method = 3;

if method == 1
    steady_state_jacobi(L, n, tolerance); %function for jacobi method
elseif method == 2
    steady_state_gauss_seidel(L, n, tolerance); % function for gauss-seidel
elseif method == 3
    steady_state_SOR(L, n, tolerance, omega); % function for SOR
end

