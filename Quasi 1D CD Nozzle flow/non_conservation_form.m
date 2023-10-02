% function to compute solution to quasi 1D subsonic-supersonic isentropic
% flow nozzle in non-conservation form
clear all
close all
clc

%function [rho,T,p,M,m,n_iter]= non_conservation_form(L,n,nt,C,tol)
L = 3;
n = 61;
throat = (n+1)/2;
x = linspace(0,3,n);
nt = 1400;
C = 0.5;
tol = 1e-4;
    % initialization
    x = linspace(0,L,n);
    dx = L/(n-1);
    A = 1 + 2.2*(x-1.5).^2;
    throat = (n+1)/2;   % assuming equal no of grid points on either side of throat
    nt_x = linspace(1,nt,nt);
    gamma = 1.4;
    error = 1;
    n_iter = 1;
    m_old = 0;

    % Initial conditions
    rho = 1 - 0.3146*x;
    T = 1 - 0.2314*x;
    v = (0.1 + 1.09.*x).*T.^(1/2);
    
    % analytical steady-state mass flow rate solution
    m_an = 0.579*ones(n);
    
    % copy of initial values
    rho_initial = rho;
    T_initial = T;
    v_initial = v;
    m_initial = rho.*v.*A;

    % loop to calculate minimum value of dt among all nodes
    for k = 1:n
        dt_x(k) = C*dx/(T(k)^(1/2)+v(k));
    end
    dt = min(dt_x);

    % time loop
    for j = 1:nt
        
        % updating old values at the start of every time loop
        rho_old = rho;
        v_old = v;
        T_old = T;

        % predictor step
        for i = 2:n-1
            drho_dx = (rho(i+1)-rho(i))/dx;
            dv_dx = (v(i+1)-v(i))/dx;
            dT_dx = (T(i+1)-T(i))/dx;
            dlnA_dx = (log(A(i+1))-log(A(i)))/dx;

            % computing the derivatives using forward difference
            % continuity equation
            drho_dt_p(i) = -rho(i)*dv_dx-rho(i)*v(i)*dlnA_dx-v(i)*drho_dx;
            % momentum equation
            dv_dt_p(i) = -v(i)*dv_dx-(1/gamma)*(dT_dx+(T(i)/rho(i))*drho_dx);
            % energy equation
            dT_dt_p(i) = -v(i)*dT_dx-(gamma-1)*T(i)*(dv_dx+v(i)*dlnA_dx);

            % predictor step update
            rho(i) = rho_old(i)+drho_dt_p(i)*dt;
            v(i) = v_old(i)+dv_dt_p(i)*dt;
            T(i) = T_old(i)+dT_dt_p(i)*dt;
        end

        % corrector step
        for i = 2:n-1
            drho_dx = (rho(i)-rho(i-1))/dx;
            dv_dx = (v(i)-v(i-1))/dx;
            dT_dx = (T(i)-T(i-1))/dx;
            dlnA_dx = (log(A(i))-log(A(i-1)))/dx;

            % computing the derivatives using rearward difference
            % continuity equation
            drho_dt_c(i) = -rho(i)*dv_dx-rho(i)*v(i)*dlnA_dx-v(i)*drho_dx;
            % momentum equation
            dv_dt_c(i) = -v(i)*dv_dx-(1/gamma)*(dT_dx+(T(i)/rho(i))*drho_dx);
            % energy equation
            dT_dt_c(i) = -v(i)*dT_dx-(gamma-1)*T(i)*(dv_dx+v(i)*dlnA_dx);
        end

        % average value of time derivatives
        drho_dt = 0.5*(drho_dt_p+drho_dt_c);
        dv_dt = 0.5*(dv_dt_p+dv_dt_c);
        dT_dt = 0.5*(dT_dt_p+dT_dt_c);
        
        % corrector step update
        for i = 2:n-1            
            rho(i) = rho_old(i)+drho_dt(i)*dt;
            v(i) = v_old(i)+dv_dt(i)*dt;
            T(i) = T_old(i)+dT_dt(i)*dt;
        end

        % Inlet boundary condition
        v(1) = 2*v(2)-v(3);

        % Outlet boundary conditions
        v(n) = 2*v(n-1)-v(n-2);
        rho(n) = 2*rho(n-1)-rho(n-2);
        T(n) = 2*T(n-1)-T(n-2);

        % Calculating pressure, mass flow rate and Mach number
        p = rho.*T;
        M = v./T.^0.5;
        m = rho.*v.*A;
        m_t(j,:) = rho.*v.*A; % capturing timewise variation of massflow at all nodes

        % Capturing flow field variables at the throat for each time step
        rho_throat(j) = rho(throat);
        T_throat(j) = T(throat);
        p_throat(j) = p(throat);
        M_throat(j) = M(throat);

        % checking convergence
        error = max(abs(m-m_old));
        if error>tol
            n_iter = n_iter + 1;            
            m_old = m;
        end

    end

##    % plotting the steady state solutions
##    figure(1);
##    sgtitle(sprintf(['Steady-state distribution of flow-field variables across nozzle'...
##        '\n-Non-conservation form, No. of time steps = %g'],nt));
##    subplot(4,1,1);
##    plot(x,M,'LineWidth',2);
##    xlabel('x/L');
##    ylabel('M');
##    title(sprintf('Mach number'));
##    subplot(4,1,2);
##    plot(x,rho,'r','LineWidth',2);
##    xlabel('x/L');
##    ylabel('\rho/\rho_o');
##    title(sprintf('Non-dimensional Density'));
##    subplot(4,1,3);
##    plot(x,T,'g','LineWidth',2);
##    xlabel('x/L');
##    ylabel('T/T_o');
##    title(sprintf('Non-dimensional Temperature'));
##    subplot(4,1,4);
##    plot(x,p,'m','LineWidth',2);
##    xlabel('x/L');
##    ylabel('p/p_o');
##    title(sprintf('Non-dimensional Pressure'));
##
##    % Timewise variarion of flow field variables at nozzle throat
##    figure(2)
##    plot(nt_x,rho_throat,'LineWidth',2);
##    hold on
##    plot(nt_x,T_throat,'LineWidth',2);
##    hold on
##    plot(nt_x,p_throat,'LineWidth',2);
##    hold on
##    plot(nt_x,M_throat,'LineWidth',2);
##    title(sprintf(['Timewise variation of non-dimensional flow-field variables'...
##        '\nat nozzle throat- Non-conservation form']));
##    xlabel('Number of time steps');
##    legend('\rho/\rho_o','T/T_o','p/p_o','M');
##
##    % Mass flow variation across nozzle at different time steps
##    figure(3)
##    plot(x,m_initial,'--',x,m_t(50,:),x,m_t(100,:),x,m_t(150,:),x,m_t(200,:),x,m_t(700,:),x,m_an,':','LineWidth',2);
##    xlabel('Non-dimensional distance through nozzle(x)');
##    ylabel('Non-dimensional mass flow (\rhoVA)/(\rho_oa_oA*)');
##    legend('0\Deltat','50\Deltat','100\Deltat','150\Deltat','200\Deltat','700\Deltat','Analytical solution');
##    title(sprintf(['Distribution of non-dimensional mass flow rate across the nozzle'...
##        '\n-Non-conservation form']));

%Calculations for graphic display start
ymax=sqrt(max(A)/pi);
ny=1001;
nx=1001;
v_cont=-0.1*ones(ny,nx);
T_cont=-0.1*ones(ny,nx);
p_cont=-0.1*ones(ny,nx);
rho_cont=-0.1*ones(ny,nx);
y=linspace(-ymax,ymax,ny);
dy=y(2)-y(1);
x_new=linspace(0,3,nx);
for i=1:nx
i_n=floor((n-1)*(i-1)/(nx-1))+1;
if nx>i
v_x(i)=v(i_n)+((v(i_n+1)-v(i_n))*((x_new(i)-x(i_n))/(x(i_n+1)-x(i_n))));
A_x(i)=A(i_n)+((A(i_n+1)-A(i_n))*((x_new(i)-x(i_n))/(x(i_n+1)-x(i_n))));
else
v_x(i)=v(i_n);
A_x(i)=A(i_n);
end
end
for i=1:nx
b=sqrt(4*A_x(i)/pi);
ny_start=round(((-b/2)+ymax)/dy)+1;
ny_end=length(y)-ny_start+1;
for j=ny_start:ny_end
v_cont(j,i)=v_x(i);
endfor
endfor
%Calculations for graphic display end
       %plot results
       figure(1)
       subplot(4,1,1)
       contourf(x_new,y,v_cont, 100, 'LineStyle', 'none')
       %colorbar
       colormap(jet)
       xlabel('Length of nozzle')
       ylabel('Breadth of nozzle')
       text_time=sprintf('Steady state time=%d s',time);
       title({'Steady-State Solution- Velocity'})
       subplot(4,1,2)
        plot(x,M,'LineWidth',2,'color','green');
        set(gca, 'color', 'k');
        xlabel('Length of nozzle')
        ylabel('Mach Number')
        subplot(4,1,3)
        plot(x,T,'LineWidth',2,'color','green');
        set(gca, 'color', 'k');
        xlabel('Length of nozzle')
        ylabel('Temp. ratio')
        subplot(4,1,4)
        plot(x,p,'LineWidth',2,'color','green');
        set(gca, 'color', 'k');
        xlabel('Length of nozzle')
        ylabel('Pressure ratio')
        
##       figure(2)
##       plot(x,v,'color','r');
##       hold on
##       plot(x,T,'color','b');
##       hold on
##       plot(x,rho,'color','g');
##       hold on
##       plot(x,p,'color','m');
##       hold off
##       xlabel('x')
##       ylabel('output ratios')
##       legend('Velocity','Temperature','Density','Pressure')
##       title_text1=sprintf('Time=%d',t);
##       title_text2=sprintf('Steady state time=%d and No. of time steps=%d',time,tn);
##       title({title_text1,title_text2});
##       axis([0 3 -0.5 2.5]);
%end
