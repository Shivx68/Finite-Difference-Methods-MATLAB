% Jacobi iterative solver

function [x, n_iter] = jacobi(T, c, tol)
    
    % Initialization   
    x = [1; 1; 1];
    x_old = x;
    error = 1;
    n_iter = 1;      
    
    % Convergence loop
    while error > tol
        x = T*x_old + c;                % Fixed point method, x_(n+1) = P*x_n + q
        error = max(abs(x - x_old));    % Convergence criterion
        x_old = x;                      % Updating old values
        n_iter = n_iter + 1;
    end
end