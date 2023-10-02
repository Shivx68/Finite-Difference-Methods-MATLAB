% function to compute solution to quasi 1D subsonic-supersonic isentropic
% flow nozzle in conservation form

 
    clear all
    close all
    clc

    L = 3;
    n = 61;
    throat = (n+1)/2;
    x = linspace(0,3,n);
    nt = 50000;
    C = 0.5;
    tol = 1e-4;
    m_analytical = 0.579*ones(n);
    % Initialization
    x = linspace(0,L,n);
    dx = L/(n-1);
    A = 1 + 2.2*(x-1.5).^2;
    throat = (n+1)/2; % assuming equal no of grid points on either side of throat
    nt_x = linspace(1,nt,nt);
    gamma = 1.4;
    error = 1;
    n_iter = 1;
    m_old = 0;

    % Initial conditions
    u2_i = 0.59;   % initial mass flow rate, at t=0;
    for i= 1:n
        if ((x(i)>=0) && (x(i)<=0.5))
            rho(i) = 1;
            T(i) = 1;
        elseif ((x(i)<1.5) && (x(i)>0.5))
            rho(i) = 1-0.366*(x(i)-0.5);
            T(i) = 1-0.167*(x(i)-0.5);
        elseif ((x(i)>=1.5) && (x(i)<=3.5))
            rho(i) = 0.634-0.3879*(x(i)-1.5);
            T(i) = 0.833-0.3507*(x(i)-1.5);
        end
        v(i) = u2_i/(rho(i)*A(i));
        u1(i) = rho(i)*A(i);
        u2(i) = rho(i)*A(i)*v(i);
        u3(i) = rho(i)*A(i)*(T(i)/(gamma-1)+(gamma/2)*v(i)^2);
    end
    
    % analytical steady-state mass flow rate solution
    m_an = 0.579*ones(n);   
    
    % copy of initial values
    u1_initial = u1;
    u2_initial = u2;
    u3_initial = u3;

    % loop to calculate minimum value of dt among all nodes
    for k = 1:n
        dt_x(k) = C*dx/(T(k)^(1/2)+v(k));
    end
    dt = min(dt_x);
    %dt = 0.0267;

    % time loop
    for j = 1:nt
        
        % updating old values at the start of every time loop
        u1_old = u1;
        u2_old = u2;
        u3_old = u3;
        
        % predictor step
        for i = 2:n-1
            
            % computing flux and source terms
            f1_i = u2(i);
            f1_ip1 = u2(i+1);
            f2_i = u2(i)^2/u1(i)+((gamma-1)/gamma)*(u3(i)-(gamma/2)*u2(i)^2/u1(i));
            f2_ip1 = u2(i+1)^2/u1(i+1)+((gamma-1)/gamma)*(u3(i+1)-(gamma/2)*u2(i+1)^2/u1(i+1));
            j2  = (1/gamma)*rho(i)*T(i)*((A(i+1) - A(i))/dx);
            f3_i = gamma*u2(i)*u3(i)/u1(i)-(gamma*(gamma-1)/2)*u2(i)^3/u1(i)^2;
            f3_ip1 = gamma*u2(i+1)*u3(i+1)/u1(i+1)-(gamma*(gamma-1)/2)*u2(i+1)^3/u1(i+1)^2;

            % computing the derivatives using forward difference
            du1_dt_p(i) = -(f1_ip1-f1_i)/dx;    % continuity eq      
            du2_dt_p(i) = j2-(f2_ip1-f2_i)/dx;  % momentum eq
            du3_dt_p(i) = -(f3_ip1-f3_i)/dx;    % energy eq

            % updating the predicted values
            u1(i) = u1(i)+du1_dt_p(i)*dt;
            u2(i) = u2(i)+du2_dt_p(i)*dt;
            u3(i) = u3(i)+du3_dt_p(i)*dt;
        end

        % Updating primitive variables after predictor step
        rho = u1./A;
        v = u2./u1;
        T = (gamma-1)*((u3./u1) - (gamma/2)*v.^2);
    
        % corrector step
        for i = 2:n-1
            
            % computing new set of flux and source terms using predicted
            % values for the dependent variables
            f1_i = u2(i);
            f1_im1 = u2(i-1);
            f2_i = u2(i)^2/u1(i)+((gamma-1)/gamma)*(u3(i)-(gamma/2)*u2(i)^2/u1(i));
            f2_im1 = u2(i-1)^2/u1(i-1)+((gamma-1)/gamma)*(u3(i-1)-(gamma/2)*u2(i-1)^2/u1(i-1));
            j2  = (1/gamma)*rho(i)*T(i)*((A(i) - A(i-1))/dx);
            f3_i = gamma*u2(i)*u3(i)/u1(i)-(gamma*(gamma-1)/2)*u2(i)^3/u1(i)^2;
            f3_im1 = gamma*u2(i-1)*u3(i-1)/u1(i-1)-(gamma*(gamma-1)/2)*u2(i-1)^3/u1(i-1)^2;

            % computing the derivatives using rearward difference
            du1_dt_c(i) = -(f1_i-f1_im1)/dx;        
            du2_dt_c(i) = j2-(f2_i-f2_im1)/dx;        
            du3_dt_c(i) = -(f3_i-f3_im1)/dx;
        end

        % computing average value of time derivatives
        du1_dt = 0.5*(du1_dt_p+du1_dt_c);
        du2_dt = 0.5*(du2_dt_p+du2_dt_c);
        du3_dt = 0.5*(du3_dt_p+du3_dt_c);

        for i = 2:n-1
            % computing corrected values
            u1(i) = u1_old(i)+du1_dt(i)*dt;
            u2(i) = u2_old(i)+du2_dt(i)*dt;
            u3(i) = u3_old(i)+du3_dt(i)*dt;
        end
        
        % applying boundary conditions
        % inlet
        u1(1) = rho(1)*A(1);
        u2(1) = 2*u2(2)-u2(3);
        v(1) = u2(1)/(rho(1)*A(1));
        u3(1) = u1(1)*(T(1)/(gamma-1)+(gamma/2)*v(1)^2);

        % outlet
        u1(n) = 2*u1(n-1)-u1(n-2);
        u2(n) = 2*u2(n-1)-u2(n-2);
        u3(n) = 2*u3(n-1)-u3(n-2);
        
        
##        % computing the values of the primitive variables
##        for i= 1:n
##            rho(i) = u1(i)/A(i);
##            v(i) = u2(i)/u1(i);
##            T(i) = (gamma-1)*(u3(i)/u1(i)-gamma/2*v(i)^2);
##        end
##
##        % Calculating pressure, mass flow rate and Mach number
##        p = rho.*T;
##        M = v./T.^0.5;
##        m = rho.*v.*A;
##        m_t(j,:) = rho.*v.*A; % capturing timewise variation of massflow at all nodes
##
##        % Capturing flow field variables at the throat for each time step
##        rho_throat(j) = rho(throat);
##        T_throat(j) = T(throat);
##        p_throat(j) = p(throat);
##        M_throat(j) = M(throat);
##
##        % checking convergence
##        error = max(abs(m-m_old));
##        if error>tol
##            n_iter = n_iter + 1;            
##            m_old = m;
          disp(j);
        end
##    end
