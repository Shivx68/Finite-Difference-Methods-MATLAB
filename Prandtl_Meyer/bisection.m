function M_act = bisection (f_act,Gamma)
    a_int = 1.1; % Left limit of the interval 
    b_int = 2.9; % Right limit of the interval
    precision = 0.0000001; % Max error
    zero_f1 = sqrt((Gamma + 1)/(Gamma - 1))*(atan(sqrt(((Gamma - 1)/(Gamma + 1))*(a_int^2 - 1)))) - (atan(sqrt((a_int^2) - 1))) - f_act; % Function used to find its zero
    zero_f2 = sqrt((Gamma + 1)/(Gamma - 1))*(atan(sqrt(((Gamma - 1)/(Gamma + 1))*(((a_int + b_int)/2)^2 - 1)))) - (atan(sqrt((((a_int + b_int)/2)^2) - 1))) - f_act;
    while ((b_int-a_int)/2 > precision)
        if (zero_f1*zero_f2 <=0)
            b_int = (a_int + b_int)/2;
        else
            a_int = (a_int + b_int)/2;
        end
        zero_f1 = sqrt((Gamma + 1)/(Gamma - 1))*(atan(sqrt(((Gamma - 1)/(Gamma + 1))*(a_int^2 - 1)))) - (atan(sqrt((a_int^2) - 1))) - f_act;
        zero_f2 = sqrt((Gamma + 1)/(Gamma - 1))*(atan(sqrt(((Gamma - 1)/(Gamma + 1))*(((a_int + b_int)/2)^2 - 1)))) - (atan(sqrt((((a_int + b_int)/2)^2) - 1))) - f_act;
    end
    M_act = (a_int + b_int)/2; % Corrected Mach number
end