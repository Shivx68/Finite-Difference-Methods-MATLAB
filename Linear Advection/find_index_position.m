% function to return the index position of a particular value
function index = find_index_position(x, xva1)
    % calculating node step size
    dx = x(2)- x(1);
    % for loop to find the closest value and return the index position
    for i = 1: length(x)
        if abs(x(i)- xva1) < dx/2     
            index = i;   
        end                                   
    end
end