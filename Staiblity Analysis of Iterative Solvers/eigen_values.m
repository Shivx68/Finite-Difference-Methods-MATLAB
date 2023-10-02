% Function to calculate eigenvalues of a 3*3 matrix

function eig_val = eigen_values(A)
    
    % The characteristic equation of a 3*3 matrix, a*y^3 - b*y^2 + c*y - d = 0
    % a = 1;
    % b = sum of the diagonal elements or tace of the matrix
    % c = sum of minors of the diagonal elements of the matrix
    % d = determinant of the matrix

    a = 1;
    b = trace(A);
    c = 0;
    % for loop to calculate the sum of minors of all diagonal elements
    for i = 1:2
        for j = 2:3
            if i ~= j
                c = c + A(i,i)*A(j,j)-A(i,j)*A(j,i);
            end
        end
    end
    d = det(A);
    % solving the cubic characteristic equation 
    eig_val = roots([a -b c -d]);
end

