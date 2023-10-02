% To solve the give system of linear equations using iterative solvers
% and to calculate eigenvalues and spectral radius of their iteration
% matrices

clear all
close all
clc

% Given system of equations
A = [5 1 2;-3 9 4;1 2 -7];
b = [10; -14; 33];

tol = 1e-4;     % convergence criterion
omega = 0.9;    % relaxation factor

% Initialization
D = zeros(3,3);
L = zeros(3,3);
U = zeros(3,3);

% Matrix Decomposition, [A] = [D] + [L] + [U]
    for j = 1:3
        for i = 1:3
            if i==j
                D(i,j) = A(i,j);
            elseif i>j
                L(i,j) = A(i,j);
            elseif i<j
                U(i,j) = A(i,j);
            end
        end
    end

% calculating the iteration matrix and constant vector for Jacobi
T_jac = -inv(D)*(L+U);
c_jac = inv(D)*b;

% calculating the iteration matrix and constant vector for Gauss-Seidel
T_gs = -inv(D+L)*U;
c_gs = inv(D+L)*b; 

% calculating the iteration matrix and constant vector for SOR
T_sor = inv((1/omega)*D+L)*((1/omega)*D-D-U);
c_sor = inv((1/omega)*D+L)*b;

% calculating eigenvalues of iteration matrix of each method
eigval_jac = eigen_values(T_jac);
eigval_gs = eigen_values(T_gs);
eigval_sor = eigen_values(T_sor);

% calculating spectral radius
rho_jac = max(abs(eigval_jac));
rho_gs = max(abs(eigval_gs));
rho_sor = max(abs(eigval_sor));

% solving the linear system using each iteration
[x_jac, n_jac] = jacobi(T_jac, c_jac, tol);
[x_gs, n_gs] = gauss_seidel(T_gs, c_gs, tol);
[x_sor, n_sor] = sor(T_sor, c_sor, tol);


% iterative solver comparison
stem(rho_jac, n_jac,'filled');
hold on
stem(rho_gs, n_gs,'filled');
hold on
stem(rho_sor, n_sor,'filled','g');
%set(gca, 'YScale', 'log') % enable if magnificationmis <1 
title('Convergence rate for Diagonal Magnification by 1');
xlabel('Spectral radius (rho)')
ylabel('No. of iterations')
legend('Jacobi method','Gauss-Seidel method','SOR method')

