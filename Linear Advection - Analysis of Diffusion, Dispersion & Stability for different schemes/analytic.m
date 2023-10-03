function [u] = analytic(x,t,v)
% function to generate the analytical function profile (@ given t)
    
    % wave travel distance
    D = v*t;
    
    % domain parameters
    n = length(x);    
    dx = x(2)-x(1);

    % intial profile parameters
    u_high = 2;
    u_low = 0;
    u_mid = 1;
    pip_start = 0.1 + D;            
    pip_mid = 0.3 + D;
    pip_end = 0.4 + D;
       
    % for loop to find the position of node for x = pip_start/pip_mid/pip_end
    for i = 1:n
        if abs(x(i)- pip_start) < dx/2     % using decision making statement
            pip_start_index = i;           % -to get the position of node
        elseif abs(x(i)-pip_end) < dx/2    % -which has the closest value to
            pip_end_index = i;             % -to the value at x at which the
        elseif abs(x(i)-pip_mid) < dx/2    % -function spikes occur.
            pip_mid_index = i;
        end                                   
    end

    % Initial Condition
    u = ones(1,n);
    for j = 1:n
      if j >= pip_start_index && j < pip_mid_index
        u(j) = u_high;
      elseif j == pip_mid_index
        u(j) = u_low;
      elseif j > pip_mid_index && j < pip_end_index
        u(j) = (u_mid - u_low)/(pip_end-pip_mid)*(x(j)-pip_mid);
      end
    end
end