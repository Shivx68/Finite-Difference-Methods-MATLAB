% %%%% ERROR in F2 term - yet to rectify
% % Prandtl_Meyer Expansion wave - 2D steady-state simulation
% 
% clear;
% close;
% clc;
% 
% % Actual Domain Dimensions
% % Height along y = 40 m
% % Length along x = 65 m
% % Expansion angle = 5.352 deg
% % Expansion point located at x = 10 m
% 
% % Defining the Computational Domain
% H_y = 40; % location of top boundary, m
% H = 1;   % height along eta axis, m
% L = 65;  % length along xi axis, m
% n_y = 41;  % no of nodes along eta axis
% eta = linspace(0, H, n_y); % divisions along y axis aka, eta axis.
% eta_stepsize = H/(n_y-1);
% theta = 5.352;
% gamma = 1.4;
% R = 287.1848; % J/kg-K
% xi = 0;
% x = xi*ones(n_y,1);
% y = eta'*H_y;
% 
% %xi_stepsize = 0.1; % only for the time being, must be determined using CFL criteria
% %n_x = round(L/xi_stepsize);
% 
% 
% % Defining the initial data line
% p = 1.01e5*ones(n_y,1); % N/m^2
% rho = 1.23*ones(n_y,1); % Kg/m^3
% T = 286.1*ones(n_y,1);  % K
% M = 2*ones(n_y,1);
% a = sqrt(gamma*R*T); % speed of sound at 286.1 K
% u = M.*a;
% v = zeros(n_y,1);
% 
% 
% Co = 0.5;
% Cy = 0.6; % Range (0.01 - 0.3)
% 
% iter = 0;
% 
% 
% 
% 
% %while (eps <= 65)
% %while i<=4
% %  i = i+1;
% 
%   % Initializing F before first marching step
% 
%   % for j = 1:n_y
%   %   F1a(j,i) = rho(j,i)*u(j,i);
%   %   F2a(j,i) = rho(j,i)*u(j,i)^2 + p(j,i);
%   %   F3a(j,i) = rho(j,i)*u(j,i)*v(j,i);
%   %   F4a(j,i) = (gamma/(gamma-1))*p(j,i)*u(j,i) + rho(j,i)*u(j,i)*(u(j,i)^2 + v(j,i)^2)/2;
%   % end
% 
%   % F1 = rho.*u;
%   % F2 = rho.*u.^2 + p;
%   % F3 = rho.*u.*v;
%   % F4 = (gamma/(gamma-1))*p.*u + rho.*u.*(u.^2 + v.^2)/2;


  
  %%%% Start of Marching loop
  for i = 29
  
  % Updating F and G for each marching step

  % for j = 1:n_y
  %   G1a(j,i) = rho(j,i)*v(j,i);
  %   G2a(j,i) = rho(j,i)*u(j,i)*v(j,i);
  %   G3a(j,i) = rho(j,i)*v(j,i)^2 + p(j,i);
  %   G4a(j,i) = gamma/(gamma-1)*p(j,i)*v(j,i) + rho(j,i)*v(j,i)*(u(j,i)^2 + v(j,i)^2)/2;
  % 
  % end
  F1(:,i) = rho(:,i).*u(:,i);
  F2(:,i) = rho(:,i).*u(:,i).^2 + p(:,i);
  F3(:,i) = rho(:,i).*u(:,i).*v(:,i);
  F4(:,i) = (gamma/(gamma-1))*p(:,i).*u(:,i) + rho(:,i).*u(:,i).*(u(:,i).^2 + v(:,i).^2)/2;

  G1(:,i) = rho(:,i).*F3(:,i)./F1(:,i);
  G2(:,i) = F3(:,i);
  G3(:,i) = rho(:,i).*(F3(:,i).^2./F1(:,i).^2) + F2(:,i) - (F1(:,i).^2./rho(:,i));
  G4(:,i) = gamma/(gamma-1)*(F2(:,i)-F1(:,i).^2./rho(:,i)).*F3(:,i)./F1(:,i) + (rho(:,i)/2).*F3(:,i)./F1(:,i).*(F1(:,i).^2./rho(:,i).^2 + F3(:,i).^2./F1(:,i).^2);

   % Adaptive marching step size
  %xi_stepsize = 0.5; % only for the time being, must be determined using CFL criteria
  %n_x = L/xi_stepzize + 1;

mu = asind(1./M(:,i));
value1 = abs(tand(theta + mu));
value2 = abs(tand(theta - mu));
value = max(value1, value2);
dy(i) = y(2,i) - y(1,i);
dx(i) = min(dy(i)./value);
  
  xi_stepsize(i) = Co*dx(i);
  %xi_stepsize(i) = 0.1;
  
  xi(i+1) = xi(i) + xi_stepsize(i);

  % Initializing h and deta_dx for each marching iteration
  %eta_a = 0;
  x(:,i+1) = xi(i+1)*ones(n_y,1);
  if (xi(i) <= 10)
    h = H_y;
    deta_dx = zeros(n_y,1);
    y_s = 0;
    y(:,i+1) = y_s + eta*h;
    
  else
    h = H_y + (xi(i)-10)*tand(theta);
    deta_dx = (1-eta)*(tand(theta)/h);
    h1 = H_y + (xi(i+1)-10)*tand(theta);
    y_s = -(xi(i+1)-10)*tand(theta);
    y(:,i+1) = y_s + eta*h;

    
    % for j = 1:n_y
    %   deta_dx_a(j) = (1-eta_a)*(tand(theta)/h);
    %   eta_a = eta_a + eta_stepsize;
    % end
  end

  

  for j = 2:n_y-1

    % Predictor Step  
    dF1_dxi_p(j,i) = deta_dx(j)*(F1(j,i) - F1(j+1,i))/eta_stepsize + (1/h)*(G1(j,i) - G1(j+1,i))/eta_stepsize;
    dF2_dxi_p(j,i) = deta_dx(j)*(F2(j,i) - F2(j+1,i))/eta_stepsize + (1/h)*(G2(j,i) - G2(j+1,i))/eta_stepsize;
    dF3_dxi_p(j,i) = deta_dx(j)*(F3(j,i) - F3(j+1,i))/eta_stepsize + (1/h)*(G3(j,i) - G3(j+1,i))/eta_stepsize;
    dF4_dxi_p(j,i) = deta_dx(j)*(F4(j,i) - F4(j+1,i))/eta_stepsize + (1/h)*(G4(j,i) - G4(j+1,i))/eta_stepsize;

  end

  % Bottom boundary

  dF1_dxi_p(1,i) = deta_dx(1)*(F1(1,i) - F1(2,i))/eta_stepsize + (1/h)*(G1(1,i) - G1(2,i))/eta_stepsize;
  dF2_dxi_p(1,i) = deta_dx(1)*(F2(1,i) - F2(2,i))/eta_stepsize + (1/h)*(G2(1,i) - G2(2,i))/eta_stepsize;
  dF3_dxi_p(1,i) = deta_dx(1)*(F3(1,i) - F3(2,i))/eta_stepsize + (1/h)*(G3(1,i) - G3(2,i))/eta_stepsize;
  dF4_dxi_p(1,i) = deta_dx(1)*(F4(1,i) - F4(2,i))/eta_stepsize + (1/h)*(G4(1,i) - G4(2,i))/eta_stepsize;

  % Top Boundary
  dF1_dxi_p(n_y,i) = deta_dx(n_y)*(F1(n_y-1,i) - F1(n_y,i))/eta_stepsize + (1/h)*(G1(n_y-1,i) - G1(n_y,i))/eta_stepsize;
  dF2_dxi_p(n_y,i) = deta_dx(n_y)*(F2(n_y-1,i) - F2(n_y,i))/eta_stepsize + (1/h)*(G2(n_y-1,i) - G2(n_y,i))/eta_stepsize;
  dF3_dxi_p(n_y,i) = deta_dx(n_y)*(F3(n_y-1,i) - F3(n_y,i))/eta_stepsize + (1/h)*(G3(n_y-1,i) - G3(n_y,i))/eta_stepsize;
  dF4_dxi_p(n_y,i) = deta_dx(n_y)*(F4(n_y-1,i) - F4(n_y,i))/eta_stepsize + (1/h)*(G4(n_y-1,i) - G4(n_y,i))/eta_stepsize;
  
  % Computing Artificial viscosity for Predictor Step %%%%%%%%%
  for j = 2:n_y-1
    SF1_p(j,i) = Cy*abs(p(j+1,i) - 2*p(j,i) + p(j-1,i))/(p(j+1,i) + 2*p(j,i) + p(j-1,i))*(F1(j+1,i) - 2*F1(j,i) + F1(j-1,i));
    SF2_p(j,i) = Cy*abs(p(j+1,i) - 2*p(j,i) + p(j-1,i))/(p(j+1,i) + 2*p(j,i) + p(j-1,i))*(F2(j+1,i) - 2*F2(j,i) + F2(j-1,i));
    SF3_p(j,i) = Cy*abs(p(j+1,i) - 2*p(j,i) + p(j-1,i))/(p(j+1,i) + 2*p(j,i) + p(j-1,i))*(F3(j+1,i) - 2*F3(j,i) + F3(j-1,i));
    SF4_p(j,i) = Cy*abs(p(j+1,i) - 2*p(j,i) + p(j-1,i))/(p(j+1,i) + 2*p(j,i) + p(j-1,i))*(F4(j+1,i) - 2*F4(j,i) + F4(j-1,i));
  end
  
  % At bundaries
  SF1_p(1,i) = 0; SF1_p(n_y,i) = 0;
  SF2_p(1,i) = 0; SF2_p(n_y,i) = 0;
  SF3_p(1,i) = 0; SF3_p(n_y,i) = 0;
  SF4_p(1,i) = 0; SF4_p(n_y,i) = 0;

  % for j = 2:n_y-1
  %   % Predictor Update
  %   F1(j,i+1) = F1(j,i) + dF1_dxi_p(j,i)*xi_stepsize + SF1(j,i);
  %   F2(j,i+1) = F2(j,i) + dF2_dxi_p(j,i)*xi_stepsize + SF2(j,i);
  %   F3(j,i+1) = F3(j,i) + dF3_dxi_p(j,i)*xi_stepsize + SF3(j,i);
  %   F4(j,i+1) = F4(j,i) + dF4_dxi_p(j,i)*xi_stepsize + SF4(j,i);
  % end

  F1(:,i+1) = F1(:,i) + dF1_dxi_p(:,i)*xi_stepsize(i) + SF1_p(:,i);
  F2(:,i+1) = F2(:,i) + dF2_dxi_p(:,i)*xi_stepsize(i) + SF2_p(:,i);
  F3(:,i+1) = F3(:,i) + dF3_dxi_p(:,i)*xi_stepsize(i) + SF3_p(:,i);
  F4(:,i+1) = F4(:,i) + dF4_dxi_p(:,i)*xi_stepsize(i) + SF4_p(:,i);

  %%%%% Here
  % for j = 2:n_y-1
  %   % Finding rho and p from F
  %   A = F3(j,i+1)^2/(2*F1(j,i+1)) - F4(j,i+1);
  %   B = gamma/(gamma-1)*F1(j,i+1)*F2(j,i+1);
  %   C = -(gamma+1)/(2*(gamma-1))*F1(j,i+1)^3;
  % 
  %   rho(j,i+1) = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  %   %p(j,i+1) = F2(i+1,j) - F1(i+1,j)^2/rho(i+1,j);
  % end

  A = F3(:,i+1).^2./(2*F1(:,i+1)) - F4(:,i+1);
  B = gamma/(gamma-1)*F1(:,i+1).*F2(:,i+1);
  C = -(gamma+1)/(2*(gamma-1))*F1(:,i+1).^3;

  rho(:,i+1) = (-B + sqrt(B.^2 - 4*A.*C))./(2*A);
  p(:,i+1) = F2(:,i+1) - F1(:,i+1).^2./rho(:,i+1);

  
  % for j = 2:n_y-1
  %   % Computing G from rho and F
  %   G1(j,i+1) = rho(j,i+1)*F3(j,i+1)/F1(j,i+1);
  %   G2(j,i+1) = F3(j,i+1);
  %   G3(j,i+1) = rho(j,i+1)*(F3(j,i+1)/F1(j,i+1))^2 + F2(j,i+1) - F1(j,i+1)^2/rho(j,i+1);
  %   G4(j,i+1) = gamma/(gamma-1)*(F2(j,i+1) - F1(j,i+1)^2/rho(j,i+1))*F3(j,i+1)/F1(j,i+1) ...
  %       + rho(j,i+1)/2*F3(j,i+1)/F1(j,i+1)*((F1(j,i+1)/rho(j,i+1))^2 + (F3(j,i+1)/F1(j,i+1))^2);
  % end

  G1(:, i+1) = rho(:,i+1).*F3(:,i+1)./F1(:,i+1);
  G2(:, i+1) = F3(:, i+1);
  G3(:, i+1) = rho(:, i+1).*(F3(:, i+1)./F1(:, i+1)).^2 + ...
      F2(:, i+1) - F1(:, i+1).^2./rho(:, i+1);
  G4(:, i+1) = gamma/(gamma-1)*(F2(:, i+1) - F1(:, i+1).^2./rho(:, i+1)).*F3(:, i+1)./F1(:, i+1) ...
        + rho(:, i+1)/2.*F3(:, i+1)./F1(:, i+1).*((F1(:, i+1)./rho(:, i+1)).^2 + (F3(:, i+1)./F1(:, i+1)).^2);

  % Corrector step
  for j = 2:n_y-1
    dF1_dxi_c(j,i) = deta_dx(j)*(F1(j-1,i+1) - F1(j,i+1))/eta_stepsize + (1/h)*(G1(j-1,i+1) - G1(j,i+1))/eta_stepsize;
    dF2_dxi_c(j,i) = deta_dx(j)*(F2(j-1,i+1) - F2(j,i+1))/eta_stepsize + (1/h)*(G2(j-1,i+1) - G2(j,i+1))/eta_stepsize;
    dF3_dxi_c(j,i) = deta_dx(j)*(F3(j-1,i+1) - F3(j,i+1))/eta_stepsize + (1/h)*(G3(j-1,i+1) - G3(j,i+1))/eta_stepsize;
    dF4_dxi_c(j,i) = deta_dx(j)*(F4(j-1,i+1) - F4(j,i+1))/eta_stepsize + (1/h)*(G4(j-1,i+1) - G4(j,i+1))/eta_stepsize;
  end

    % Bottom boundary

  dF1_dxi_c(1,i) = deta_dx(1)*(F1(1,i+1) - F1(2,i+1))/eta_stepsize + (1/h)*(G1(1,i+1) - G1(2,i+1))/eta_stepsize;
  dF2_dxi_c(1,i) = deta_dx(1)*(F2(1,i+1) - F2(2,i+1))/eta_stepsize + (1/h)*(G2(1,i+1) - G2(2,i+1))/eta_stepsize;
  dF3_dxi_c(1,i) = deta_dx(1)*(F3(1,i+1) - F3(2,i+1))/eta_stepsize + (1/h)*(G3(1,i+1) - G3(2,i+1))/eta_stepsize;
  dF4_dxi_c(1,i) = deta_dx(1)*(F4(1,i+1) - F4(2,i+1))/eta_stepsize + (1/h)*(G4(1,i+1) - G4(2,i+1))/eta_stepsize;

  % Top Boundary
  dF1_dxi_c(n_y,i) = deta_dx(n_y)*(F1(n_y-1,i+1) - F1(n_y,i+1))/eta_stepsize + (1/h)*(G1(n_y-1,i+1) - G1(n_y,i+1))/eta_stepsize;
  dF2_dxi_c(n_y,i) = deta_dx(n_y)*(F2(n_y-1,i+1) - F2(n_y,i+1))/eta_stepsize + (1/h)*(G2(n_y-1,i+1) - G2(n_y,i+1))/eta_stepsize;
  dF3_dxi_c(n_y,i) = deta_dx(n_y)*(F3(n_y-1,i+1) - F3(n_y,i+1))/eta_stepsize + (1/h)*(G3(n_y-1,i+1) - G3(n_y,i+1))/eta_stepsize;
  dF4_dxi_c(n_y,i) = deta_dx(n_y)*(F4(n_y-1,i+1) - F4(n_y,i+1))/eta_stepsize + (1/h)*(G4(n_y-1,i+1) - G4(n_y,i+1))/eta_stepsize;

  % 
  % for j = 2:n_y-1
  %   % Finding the average derivatives
  %   dF1_deps_avg(i,j) = 1/2*(dF1_deps_p(i,j)+dF1_deps_c(i,j));
  %   dF2_deps_avg(i,j) = 1/2*(dF2_deps_p(i,j)+dF2_deps_c(i,j));
  %   dF3_deps_avg(i,j) = 1/2*(dF3_deps_p(i,j)+dF3_deps_c(i,j));
  %   dF4_deps_avg(i,j) = 1/2*(dF4_deps_p(i,j)+dF4_deps_c(i,j));
  % end

   dF1_dxi_avg(:,i) = 0.5*(dF1_dxi_p(:,i) + dF1_dxi_c(:,i));
   dF2_dxi_avg(:,i) = 0.5*(dF2_dxi_p(:,i) + dF2_dxi_c(:,i));
   dF3_dxi_avg(:,i) = 0.5*(dF3_dxi_p(:,i) + dF3_dxi_c(:,i));
   dF4_dxi_avg(:,i) = 0.5*(dF4_dxi_p(:,i) + dF4_dxi_c(:,i));

  % 
  % Computing Artificial Viscosity for the corrector step
  for j = 2:n_y-1 
    SF1_c(j,i) = Cy*abs(p(j+1,i+1) - 2*p(j,i+1) + p(j-1,i+1))/(p(j+1,i+1) + 2*p(j,i+1) + p(j-1,i+1))*(F1(j+1,i+1) - 2*F1(j,i+1) + F1(j-1,i+1));
    SF2_c(j,i) = Cy*abs(p(j+1,i+1) - 2*p(j,i+1) + p(j-1,i+1))/(p(j+1,i+1) + 2*p(j,i+1) + p(j-1,i+1))*(F2(j+1,i+1) - 2*F2(j,i+1) + F2(j-1,i+1));
    SF3_c(j,i) = Cy*abs(p(j+1,i+1) - 2*p(j,i+1) + p(j-1,i+1))/(p(j+1,i+1) + 2*p(j,i+1) + p(j-1,i+1))*(F3(j+1,i+1) - 2*F3(j,i+1) + F3(j-1,i+1));
    SF4_c(j,i) = Cy*abs(p(j+1,i+1) - 2*p(j,i+1) + p(j-1,i+1))/(p(j+1,i+1) + 2*p(j,i+1) + p(j-1,i+1))*(F4(j+1,i+1) - 2*F4(j,i+1) + F4(j-1,i+1));
  end

  % At bundaries
  SF1_c(1,i) = 0; SF1_c(n_y,i) = 0;
  SF2_c(1,i) = 0; SF2_c(n_y,i) = 0;
  SF3_c(1,i) = 0; SF3_c(n_y,i) = 0;
  SF4_c(1,i) = 0; SF4_c(n_y,i) = 0;

  % 
  % for j = 2:n_y-1  
  %   % Final Update
  %   F1(i+1,j) = F1(i,j) + dF1_deps_avg(i,j)*deps + SF1(i+1,j);
  %   F2(i+1,j) = F2(i,j) + dF2_deps_avg(i,j)*deps + SF2(i+1,j);
  %   F3(i+1,j) = F3(i,j) + dF3_deps_avg(i,j)*deps + SF3(i+1,j);
  %   F4(i+1,j) = F4(i,j) + dF4_deps_avg(i,j)*deps + SF4(i+1,j);
  % end

  F1(:,i+1) = F1(:,i) + dF1_dxi_avg(:,i)*xi_stepsize(i) + SF1_c(:,i);
  F2(:,i+1) = F2(:,i) + dF2_dxi_avg(:,i)*xi_stepsize(i) + SF2_c(:,i);
  F3(:,i+1) = F3(:,i) + dF3_dxi_avg(:,i)*xi_stepsize(i) + SF3_c(:,i);
  F4(:,i+1) = F4(:,i) + dF4_dxi_avg(:,i)*xi_stepsize(i) + SF4_c(:,i);

  
  % 
  % % Final Update of the primitive variables
  % for j = 2:n_y-1
  %   % Finding corrected rho from corrected F
  %   A = F3(i+1,j)^2/(2*F1(i+1,j)) - F4(i+1,j);
  %   B = gamma/(gamma-1)*F1(i+1,j)*F2(i+1,j);
  %   C = -(gamma+1)/(2*(gamma-1))*F1(i+1,j)^3;
  % 
  %   rho(i+1,j) = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  % end

  A = F3(:,i+1).^2./(2*F1(:,i+1)) - F4(:,i+1);
  B = gamma/(gamma-1)*F1(:,i+1).*F2(:,i+1);
  C = -(gamma+1)/(2*(gamma-1))*F1(:,i+1).^3;

  rho(:,i+1) = (-B + sqrt(B.^2 - 4*A.*C))./(2*A);
  p(:,i+1) = F2(:,i+1) - F1(:,i+1).^2./rho(:,i+1);
  u(:,i+1) = F1(:,i+1)./rho(:,i+1);
  v(:,i+1) = F3(:,i+1)./F1(:,i+1);
  T(:,i+1) = p(:,i+1)./(rho(:,i+1)*R);
  a(:,i+1) = sqrt(gamma*R*T(:,i+1));
  M(:,i+1) = (sqrt(u(:,i+1).^2 + v(:,i+1).^2))./a(:,i+1);

  % Applying Slip BC at the wall - Using Prandtl-Meyer function to turn the
  % flow
  rho_cal1 = rho(1,i+1);
  u_cal1 = u(1,i+1);
  v_cal1 = v(1,i+1);
  p_cal1 = p(1,i+1); 
  T_cal1 = T(1,i+1); 
  a_cal1 = a(1,i+1);
  M_cal1 = M(1,i+1);

   
  % Computing the calculated Prandtl-Meyer function
  f_cal1 = prandtl_meyer_function(M_cal1,gamma);

  % Computing the Prandtl-Meyer rotation angle
  if (xi(i) <= 10)
    phi1 = atand(v_cal1/u_cal1);
  else
    psi1 = atand(abs(v_cal1)/u_cal1);
    phi1 = theta - psi1;
  end

  % 
  % Computing the actual Prandtl-Meyer function
  f_act1 = f_cal1 + phi1;

  M_act1 = fsolve(@(M) F (M,f_act1,gamma), M_cal1);
  % 

  % % Computing the M_act using Newton-Raphson method
  % tol = 1e-4;
  % error = 1;
  % RF = 0.1;
  % dM = 1e-4;
  % %M_act = 1;  % guess value
  % while error > tol
  %   M_act1a = M_cal1 - RF*(F(M_cal1, f_act1, gamma)/F_prime(f_act1, M_cal1, dM));
  %   error = abs(F(M_act1a, f_act1, gamma));
  %   M_cal1 = M_act1a;
  % end
  % 

  % Computing the actual values of primitive variables on the boundary
  p_act1 = p_cal1*((1 + ((gamma-1)/2)*M_cal1^2)/(1 + ((gamma-1)/2)*M_act1^2))^(gamma/(gamma-1));
  T_act1 = T_cal1*((1 + ((gamma-1)/2)*M_cal1^2)/(1 + ((gamma-1)/2)*M_act1^2));
  rho_act1 = p_act1/(R*T_act1);
  % 
  % Updating the boundary node with actual boundary values
  p(1,i+1) = p_act1;
  T(1,i+1) = T_act1;
  rho(1,i+1) = rho_act1;
  u(1,i+1) = u_cal1;
  M(1,i+1) = M_act1;
  a(1,i+1) = (sqrt(gamma*R*T_act1));

  if (xi(i) <= 10)
       v(1,i+1) = 0;
  else
      v(1,i+1) = - u_cal1*tand(theta);
  end
  %

  % % Top boundary
  rho_calny = rho(n_y,i+1);
  u_calny = u(n_y,i+1);
  v_calny = v(n_y,i+1);
  p_calny = p(n_y,i+1); 
  T_calny = T(n_y,i+1); 
  a_calny = a(n_y,i+1);
  M_calny = M(n_y,i+1);

   
  % Computing the calculated Prandtl-Meyer function
  f_calny = prandtl_meyer_function(M_calny,gamma);

  % Computing the Prandtl-Meyer rotation angle
  if (xi <= 10)
    phiny = atand(v_calny/u_calny);
  else
    psiny = atand(abs(v_calny)/u_calny);
    phiny = theta - psiny;
  end

  % 
  % Computing the actual Prandtl-Meyer function
  f_actny = f_calny + phiny;

  M_actny = fsolve(@(M) F (M,f_actny,gamma), M_calny);
  
  % Computing the actual values of primitive variables on the boundary
  p_actny = p_calny*((1 + ((gamma-1)/2)*M_calny^2)/(1 + ((gamma-1)/2)*M_actny^2))^(gamma/(gamma-1));
  T_actny = T_calny*((1 + ((gamma-1)/2)*M_calny^2)/(1 + ((gamma-1)/2)*M_actny^2));
  rho_actny = p_actny/(R*T_actny);

  % Updating the boundary node with actual boundary values
  p(n_y,i+1) = p_actny;
  T(n_y,i+1) = T_actny;
  rho(n_y,i+1) = rho_actny;
  u(n_y,i+1) = u_calny;
  M(n_y,i+1) = M_actny;
  a(n_y,i+1) = (sqrt(gamma*R*T_actny));
  if (xi(i) <= 10)
       v(n_y,i+1) = 0;
  else
      v(n_y,i+1) = - u_calny*tand(theta);
  end

  
  end
contourf(x,y,M);
  colorbar;
  %pause(0.1);



















