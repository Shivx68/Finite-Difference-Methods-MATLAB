clear all
close all
clc
%inputs
n=31;
x=linspace(0,3,n);
dx=x(2)-x(1);
gamma=1.4;
%calculate initial profile
rho=1-0.3146*x;
T=1-0.2314*x;
v=(0.1+1.09*x).*T.^0.5;
p=rho.*T;
%area
A=1+2.2*(x-1.5).^2;
%Time loop
nt=5000;
t=0;%to find the net time
%to find the time at which system stabilizes
tol=1e-4;
time=0;%to find time st which system stabilizes
tn=0;%to find time step number st which system stabilizes
C=0.5;%couriant Number
%outer time loop
for k=1:nt
%dt change as stated in anderson textbook
delta_t=C.*(dx./((T.^0.5)+v));
dt=min(delta_t);
%calculating time
t=t+dt;
%storing old solution vector values
    rho_old=rho;
    v_old=v;
    T_old=T;
    p_old=p;
        %predictor step
        for j=2:n-1
        dvdx=(v(j+1)-v(j))/dx;
        drhodx=(rho(j+1)-rho(j))/dx;
        dTdx=(T(j+1)-T(j))/dx;
        d_logA_dx=(log(A(j+1))-log(A(j)))/dx;
          %Continuity Equation
          drhodt_pred(j)=-rho(j)*dvdx-rho(j)*v(j)*d_logA_dx-v(j)*drhodx;
          %momentum equation
          dvdt_pred(j)=-v(j)*dvdx-(1/gamma)*(dTdx+(T(j)/rho(j))*drhodx);
          %energy equation
          dTdt_pred(j)=-v(j)*dTdx-(gamma-1)*T(j)*(dvdx+v(j)*d_logA_dx);
          %solution update
          rho(j)=rho(j)+drhodt_pred(j)*dt;
          v(j)=v(j)+dvdt_pred(j)*dt;
          T(j)=T(j)+dTdt_pred(j)*dt;
           endfor
         %Corrextor method
        for j=2:n-1
        dvdx=(v(j)-v(j-1))/dx;
        drhodx=(rho(j)-rho(j-1))/dx;
        dTdx=(T(j)-T(j-1))/dx;
        d_logA_dx=(log(A(j))-log(A(j-1)))/dx;
          %Continuity Equation
          drhodt_c(j)=-rho(j)*dvdx-rho(j)*v(j)*d_logA_dx-v(j)*drhodx;
          %momentum equation
          dvdt_c(j)=-v(j)*dvdx-(1/gamma)*(dTdx+(T(j)/rho(j))*drhodx);
          %energy equation
          dTdt_c(j)=-v(j)*dTdx-(gamma-1)*T(j)*(dvdx+v(j)*d_logA_dx);
         endfor
         % compute average derivatives
         drhodt=0.5*(drhodt_c+drhodt_pred);
         dvdt=0.5*(dvdt_c+dvdt_pred);
         dTdt=0.5*(dTdt_c+dTdt_pred);
         %Final Solution
         for i=2:n-1
         rho(i) =rho_old(i)+drhodt(i)*dt;
         v(i) =v_old(i)+dvdt(i)*dt;
         T(i) =T_old(i)+dTdt(i)*dt;
       endfor
       %apply boundary condition
       %inlet
       v(1)=2*v(2)-v(3);
       %outlet
       v(n)=2*v(n-1)-v(n-2);
      rho(n)=2*rho(n-1)-rho(n-2);
       T(n)=2*T(n-1)-T(n-2);
      %pressure calculation
       p=rho.*T;
       %Steady state coding to find time at which solution stabilizes
       a(k,1:n)=v;
       b(k,1:n)=rho;
       c(k,1:n)=T;
       d(k,1:n)=p;
       k_s=10;%time step difference to find error
       if k&gt;k_s
       error_v=max(abs(a(k,1:n)-a(k-k_s,1:n)));
       error_rho=max(abs(b(k,1:n)-b(k-k_s,1:n)));
       error_T=max(abs(c(k,1:n)-c(k-k_s,1:n)));
       error_p=max(abs(d(k,1:n)-d(k-k_s,1:n)));
       if time==0
       if tol&gt;=error_v
         if tol&gt;=error_rho
           if tol&gt;=error_T
             if tol&gt;=error_p
              time=t;
              tn=k;
              end
            end
         end
       end=20
     end
     end
     %steady state coding end
endfor
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
if nx&gt;i
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
       contourf(x_new,y,v_cont, 100)
       colorbar
       colormap(jet)
       xlabel('Length of the Convergence Divergence Nozzle')
       ylabel('Breadth')
       text_time=sprintf('Steady state time=%d s',time);
       title({'Velocity profile contour in the approximate shape of nozzle'=
,text_time})
       figure(2)
       plot(x,v,'color','r');
       hold on
       plot(x,T,'color','b');
       hold on
       plot(x,rho,'color','g');
       hold on
       plot(x,p,'color','m');
       hold off
       xlabel('x')
       ylabel('output ratios')
       legend('Velocity','Temperature','Density','Pressure')
       title_text1=sprintf('Time=%d',t);
       title_text2=sprintf('Steady state time=%d and No. of time steps=
=%d',time,tn);
       title({title_text1,title_text2});
       axis([0 3 -0.5 2.5]);