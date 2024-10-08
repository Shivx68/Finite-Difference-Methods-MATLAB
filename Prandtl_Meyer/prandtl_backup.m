%%%% ERROR in F2 term - yet to rectify
% Prandtl_Meyer Expansion wave - 2D steady-state simulation

clear all;
close all;
clc;

% Actual Domain Dimensions
% Height along y = 40 m
% Length along x = 65 m
% Expansion angle = 5.352 deg
% Expansion point located at x = 10 m

% Defining the Computational Domain
H = 1;   % height along eta axis, m
L = 65;  % length along xi axis, m
n_y = 41;  % no of nodes along eta axis
eta = linspace(0, H, n_y); % divisions along y axis aka, eta axis.
eta_stepsize = H/(n_y-1);
theta = 5.352;
gamma = 1.4;
R = 287.1848; % J/kg-K
xi = 0;



% Defining the initial data line
p = 1.01e5*ones(n_y,1); % N/m^2
rho = 1.23*ones(n_y,1); % Kg/m^3
T = 286.1*ones(n_y,1);  % K
M = 2*ones(n_y,1);
a = sqrt(gamma*R*T); % speed of sound at 286.1 K
u = M.*a;
v = zeros(n_y,1);


Co = 0.5;
Cy = 0.6; % Range (0.01 - 0.3)

iter = 0;




%while (eps <= 65)
%while i<=4
%  i = i+1;

  % Initializing F before first marching step

  % for j = 1:n_y
  %   F1a(j,i) = rho(j,i)*u(j,i);
  %   F2a(j,i) = rho(j,i)*u(j,i)^2 + p(j,i);
  %   F3a(j,i) = rho(j,i)*u(j,i)*v(j,i);
  %   F4a(j,i) = (gamma/(gamma-1))*p(j,i)*u(j,i) + rho(j,i)*u(j,i)*(u(j,i)^2 + v(j,i)^2)/2;
  % end

  F1 = rho.*u;
  F2 = rho.*u.^2 + p;
  F3 = rho.*u.*v;
  F4 = (gamma/(gamma-1))*p.*u + rho.*u.*(u.^2 + v.^2)/2;


  
  %%%% Start of Marching loop
  for i = 1
  
  % Updating G for each marching step

  % for j = 1:n_y
  %   G1a(j,i) = rho(j,i)*v(j,i);
  %   G2a(j,i) = rho(j,i)*u(j,i)*v(j,i);
  %   G3a(j,i) = rho(j,i)*v(j,i)^2 + p(j,i);
  %   G4a(j,i) = gamma/(gamma-1)*p(j,i)*v(j,i) + rho(j,i)*v(j,i)*(u(j,i)^2 + v(j,i)^2)/2;
  % 
  % end

  G1 = rho.*F3./F1;
  G2 = F3;
  G3 = rho.*(F3.^2./F1.^2) + F2 - (F1.^2./rho);
  G4 = gamma/(gamma-1)*(F2-F1.^2./rho).*F3./F1 + (rho/2).*F3./F1.*(F1.^2./rho.^2 + F3.^2./F1.^2);

  % Adaptive marching step size
  xi_stepsize = 0.1; % only for the time being, must be determined using CFL criteria
  %n_x = L/xi_stepzize + 1;

  % for j = 1:n_y-1
  %   mu = asind(1/M(i,j));
  %   value(1) = abs(tan(theta + mu));
  %   value(2) = abs(tan(theta - mu));
  %   deps_array(i,j) = Co*deta/max(value);
  % end
  % deps = min(deps_array(i,:));
  % eps = eps + deps;
  
  xi = xi + xi_stepsize;
  
  % Initializing h and deta_dx for each marching iteration
  eta_a = 0;
  if (xi <= 10)
    h = 40;
    deta_dx = zeros(n_y,1);
  else
    h = 40 + (xi-10)*tand(theta);
    deta_dx = (1-eta)*(tand(theta)/h);
    
    % for j = 1:n_y
    %   deta_dx_a(j) = (1-eta_a)*(tand(theta)/h);
    %   eta_a = eta_a + eta_stepsize;
    % end
  end

  for j = 1:n_y-1

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

  F1(:,i+1) = F1(:,i) + dF1_dxi_p(:,i)*xi_stepsize + SF1_p(:,i);
  F2(:,i+1) = F2(:,i) + dF2_dxi_p(:,i)*xi_stepsize + SF2_p(:,i);
  F3(:,i+1) = F3(:,i) + dF3_dxi_p(:,i)*xi_stepsize + SF3_p(:,i);
  F4(:,i+1) = F4(:,i) + dF4_dxi_p(:,i)*xi_stepsize + SF4_p(:,i);

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

  F1(:,i+1) = F1(:,i) + dF1_dxi_avg(:,i)*eta_stepsize + SF1_c(:,i);
  F2(:,i+1) = F2(:,i) + dF2_dxi_avg(:,i)*eta_stepsize + SF2_c(:,i);
  F3(:,i+1) = F3(:,i) + dF3_dxi_avg(:,i)*eta_stepsize + SF3_c(:,i);
  F4(:,i+1) = F4(:,i) + dF4_dxi_avg(:,i)*eta_stepsize + SF4_c(:,i);

  
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
  p_cal1 = p(1,i+1); %%%%%%%%%%%
  T_cal1 = T(1,i+1); %%%%%%%% define R
  a_cal1 = a(1,i+1);
  M_cal1 = M(1,i+1);

  % for j = 2:n_y-1
  %   % Finding corrected u, v, p, T, M from corrected rho and F
  %   ua(i+1,j) = F1(i+1,j)/rho(i+1,j);
  %   va(i+1,j) = F3(i+1,j)/F1(i+1,j);
  %   pa(i+1,j) = F2(i+1,j) - u(i+1,j)*F1(i+1,j);
  %   Ta(i+1,j) = p(i+1,j)/(rho(i+1,j)*R); %%%%%%%% define R
  %   aa(i+1,j) = sqrt(gamma*R*T(i+1,j));
  %   Ma(i+1,j) = (sqrt(v(i+1,j)^2 + u(i+1, j)^2))/a(i+1,j);
  % end
  % 
  % % Applying the boundary conditions
  % % Bottom boundary
  % % Predictor step
  % dF1_deps_p(i,1) = deta_dx(1)*(F1(i,1) - F1(i,2))/deta + (1/h)*(G1(i,1) - G1(i,2))/deta;
  % dF2_deps_p(i,1) = deta_dx(1)*(F2(i,1) - F2(i,2))/deta + (1/h)*(G2(i,1) - G2(i,2))/deta;
  % dF3_deps_p(i,1) = deta_dx(1)*(F3(i,1) - F3(i,2))/deta + (1/h)*(G3(i,1) - G3(i,2))/deta;
  % dF4_deps_p(i,1) = deta_dx(1)*(F4(i,1) - F4(i,2))/deta + (1/h)*(G4(i,1) - G4(i,2))/deta;
  % 
  % % Predictor Update
  % F1(i+1,1) = F1(i,1) + dF1_deps_p(i,1)*deps;
  % F2(i+1,1) = F2(i,1) + dF2_deps_p(i,1)*deps;
  % F3(i+1,1) = F3(i,1) + dF3_deps_p(i,1)*deps;
  % F4(i+1,1) = F4(i,1) + dF4_deps_p(i,1)*deps;
  % 
  % % Computing rho
  % A = F3(i+1,1)^2/(2*F1(i+1,1)) - F4(i+1,1);
  % B = gamma/(gamma-1)*F1(i+1,1)*F2(i+1,1);
  % C = -(gamma+1)/(2*(gamma-1))*F1(i+1,1)^3;
  % 
  % rho(i+1,1) = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  % 
  % % Computing G
  % G1(i+1,1) = rho(i+1,1)*(F3(i+1,1)/F1(i+1,1));
  % G2(i+1,1) = F3(i+1,1);
  % G3(i+1,1) = rho(i+1,1)*(F3(i+1,1)/F1(i+1,1))^2 + F2(i+1,1) - F1(i+1,1)^2/rho(i+1,1);
  % G4(i+1,1) = gamma/(gamma-1)*(F2(i+1,1) - F1(i+1,1)^2/rho(i+1,1))*F3(i+1,1)/F1(i+1,1) + rho(i+1,1)/2*F3(i+1,1)/F1(i+1,1)*((F1(i+1,1)/rho(i+1,1))^2 + (F3(i+1,1)/F1(i+1,1))^2);
  % 
  % % Corrector step
  % dF1_deps_c(i,1) = deta_dx(1)*(F1(i+1,1) - F1(i+1,2))/deta + (1/h)*(G1(i+1,1) - G1(i+1,2))/deta;
  % dF2_deps_c(i,1) = deta_dx(1)*(F2(i+1,1) - F2(i+1,2))/deta + (1/h)*(G2(i+1,1) - G2(i+1,2))/deta;
  % dF3_deps_c(i,1) = deta_dx(1)*(F3(i+1,1) - F3(i+1,2))/deta + (1/h)*(G3(i+1,1) - G3(i+1,2))/deta;
  % dF4_deps_c(i,1) = deta_dx(1)*(F4(i+1,1) - F4(i+1,2))/deta + (1/h)*(G4(i+1,1) - G4(i+1,2))/deta;
  % 
  % % Computing average
  % dF1_deps_avg(i,1) = 1/2*(dF1_deps_p(i,1)+dF1_deps_c(i,1));
  % dF2_deps_avg(i,1) = 1/2*(dF2_deps_p(i,1)+dF2_deps_c(i,1));
  % dF3_deps_avg(i,1) = 1/2*(dF3_deps_p(i,1)+dF3_deps_c(i,1));
  % dF4_deps_avg(i,1) = 1/2*(dF4_deps_p(i,1)+dF4_deps_c(i,1));
  % 
  % % Final boundary node update
  % F1(i+1,1) = F1(i,1) + dF1_deps_avg(i,1)*deps;
  % F2(i+1,1) = F2(i,1) + dF2_deps_avg(i,1)*deps;
  % F3(i+1,1) = F3(i,1) + dF3_deps_avg(i,1)*deps;
  % F4(i+1,1) = F4(i,1) + dF4_deps_avg(i,1)*deps;
  % 
  % % Computing corrected rho at the boundary
  % A = F3(i+1,1)^2/(2*F1(i+1,1)) - F4(i+1,1);
  % B = gamma/(gamma-1)*F1(i+1,1)*F2(i+1,1);
  % C = -(gamma+1)/(2*(gamma-1))*F1(i+1,1)^3;
  % 
  % rho_cal = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  % 
  % % Computing primitive boundary values
  % u_cal = F1(i+1,1)/rho_cal;
  % v_cal = F3(i+1,1)/F1(i+1,1);
  % p_cal = F2(i+1,1) - u_cal*F1(i+1,1); %%%%%%%%%%%
  % T_cal = p_cal/(rho_cal*R); %%%%%%%% define R
  % a_cal = sqrt(gamma*R*T_cal);
  % M_cal = sqrt( u_cal^2 + v_cal^2)/a_cal;
  % 
  % % Computing the calculated Prandtl-Meyer function
  % f_cal = prandtl_meyer_function(M_cal);
  % 
  % % Computing the Prandtl-Meyer rotation angle
  % if (eps <= 10)
  %   phi = atand(v_cal/u_cal);
  % else
  %   psi = atand(abs(v_cal)/u_cal);
  %   phi = theta - psi;
  % end
  % 
  % % Computing the actual Prandtl-Meyer function
  % f_act = f_cal + phi;
  % 
  % % Computing the M_act using Newton-Raphson method
  % tol = 1e-4;
  % error = 1;
  % RF = 0.1;
  % dM = 1e-4;
  % M_act = 1;  % guess value
  % while error > tol
  %   M_act = M_act - RF*(F(f_act, M_act)/F_prime(f_act, M_act, dM));
  %   error = abs(F(f_act, M_act));
  % end
  % 
  % % Computing the actual values of primitive variables on the boundary
  % p_act = p_cal*((1 + ((gamma-1)/2)*M_cal^2)/(1 + ((gamma-1)/2)*M_act^2))^(gamma/(gamma-1));
  % T_act = T_cal*((1 + ((gamma-1)/2)*M_cal^2)/(1 + ((gamma-1)/2)*M_act^2));
  % rho_act = p_act/(R*T_act);
  % 
  % % Updating the boundary node with actual boundary values
  % p(i+1,1) = p_act;
  % T(i+1,1) = T_act;
  % rho(i+1,1) = rho_act;
  % u(i+1,1) = u_cal;
  % v(i+1,1) = - u_cal*tand(theta);
  % M(i+1,1) = M_act;
  % a(i+1,1) = (sqrt(gamma*R*T(i+1,1)));
  % 
  % % Top boundary
  % % Predictor step
  % dF1_deps_p(i,n_y) = deta_dx(n_y)*(F1(i,n_y-1) - F1(i,n_y))/deta + (1/h)*(G1(i,n_y-1) - G1(i,n_y))/deta;
  % dF2_deps_p(i,n_y) = deta_dx(n_y)*(F2(i,n_y-1) - F2(i,n_y))/deta + (1/h)*(G2(i,n_y-1) - G2(i,n_y))/deta;
  % dF3_deps_p(i,n_y) = deta_dx(n_y)*(F3(i,n_y-1) - F3(i,n_y))/deta + (1/h)*(G3(i,n_y-1) - G3(i,n_y))/deta;
  % dF4_deps_p(i,n_y) = deta_dx(n_y)*(F4(i,n_y-1) - F4(i,n_y))/deta + (1/h)*(G4(i,n_y-1) - G4(i,n_y))/deta;
  % 
  % % Predictor Update
  % F1(i+1,n_y) = F1(i,n_y) + dF1_deps_p(i,n_y)*deps;
  % F2(i+1,n_y) = F2(i,n_y) + dF2_deps_p(i,n_y)*deps;
  % F3(i+1,n_y) = F3(i,n_y) + dF3_deps_p(i,n_y)*deps;
  % F4(i+1,n_y) = F4(i,n_y) + dF4_deps_p(i,n_y)*deps;
  % 
  % % Computing rho
  % A = F3(i+1,n_y)^2/(2*F1(i+1,n_y)) - F4(i+1,n_y);
  % B = gamma/(gamma-1)*F1(i+1,n_y)*F2(i+1,n_y);
  % C = -(gamma+1)/(2*(gamma-1))*F1(i+1,n_y)^3;
  % 
  % rho(i+1,n_y) = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  % 
  % % Computing G
  % G1(i+1,n_y) = rho(i+1,n_y)*(F3(i+1,n_y)/F1(i+1,n_y));
  % G2(i+1,n_y) = F3(i+1,n_y);
  % G3(i+1,n_y) = rho(i+1,n_y)*(F3(i+1,n_y)/F1(i+1,n_y))^2 + F2(i+1,n_y) - F1(i+1,n_y)^2/rho(i+1,n_y);
  % G4(i+1,n_y) = gamma/(gamma-1)*(F2(i+1,n_y) - F1(i+1,n_y)^2/rho(i+1,n_y))*F3(i+1,n_y)/F1(i+1,n_y) + rho(i+1,n_y)/2*F3(i+1,n_y)/F1(i+1,n_y)*((F1(i+1,n_y)/rho(i+1,n_y))^2 + (F3(i+1,n_y)/F1(i+1,n_y))^2);
  % 
  % % Corrector step
  % dF1_deps_c(i,n_y) = deta_dx(n_y)*(F1(i+1,n_y-1) - F1(i+1,n_y))/deta + (1/h)*(G1(i+1,n_y-1) - G1(i+1,n_y))/deta;
  % dF2_deps_c(i,n_y) = deta_dx(n_y)*(F2(i+1,n_y-1) - F2(i+1,n_y))/deta + (1/h)*(G2(i+1,n_y-1) - G2(i+1,n_y))/deta;
  % dF3_deps_c(i,n_y) = deta_dx(n_y)*(F3(i+1,n_y-1) - F3(i+1,n_y))/deta + (1/h)*(G3(i+1,n_y-1) - G3(i+1,n_y))/deta;
  % dF4_deps_c(i,n_y) = deta_dx(n_y)*(F4(i+1,n_y-1) - F4(i+1,n_y))/deta + (1/h)*(G4(i+1,n_y-1) - G4(i+1,n_y))/deta;
  % 
  % % Computing average
  % dF1_deps_avg(i,n_y) = 1/2*(dF1_deps_p(i,n_y)+dF1_deps_c(i,n_y));
  % dF2_deps_avg(i,n_y) = 1/2*(dF2_deps_p(i,n_y)+dF2_deps_c(i,n_y));
  % dF3_deps_avg(i,n_y) = 1/2*(dF3_deps_p(i,n_y)+dF3_deps_c(i,n_y));
  % dF4_deps_avg(i,n_y) = 1/2*(dF4_deps_p(i,n_y)+dF4_deps_c(i,n_y));
  % 
  % % Final boundary node update
  % F1(i+1,n_y) = F1(i,n_y) + dF1_deps_avg(i,n_y)*deps;
  % F2(i+1,n_y) = F2(i,n_y) + dF2_deps_avg(i,n_y)*deps;
  % F3(i+1,n_y) = F3(i,n_y) + dF3_deps_avg(i,n_y)*deps;
  % F4(i+1,n_y) = F4(i,n_y) + dF4_deps_avg(i,n_y)*deps;
  % 
  % % Computing corrected rho at the boundary
  % A = F3(i+1,n_y)^2/(2*F1(i+1,n_y)) - F4(i+1,n_y);
  % B = gamma/(gamma-1)*F1(i+1,n_y)*F2(i+1,n_y);
  % C = -(gamma+1)/(2*(gamma-1))*F1(i+1,n_y)^3;
  % 
  % rho_cal = (-B + sqrt(B^2 - 4*A*C))/(2*A);
  % 
  % % Computing primitive boundary values
  % u_cal = F1(i+1,n_y)/rho_cal;
  % v_cal = F3(i+1,n_y)/F1(i+1,n_y);
  % p_cal = F2(i+1,n_y) - u_cal*F1(i+1,n_y); %%%%%%%%%%%
  % T_cal = p_cal/(rho_cal*R); %%%%%%%% define R
  % a_cal = sqrt(gamma*R*T_cal);
  % M_cal = sqrt( u_cal^2 + v_cal^2)/a_cal;
  % 
  % % Computing the calculated Prandtl-Meyer function
  % f_cal = prandtl_meyer_function(M_cal);
  % 
  % % Computing the Prandtl-Meyer rotation angle
  % phi = atand(v_cal/u_cal);
  % 
  % % Computing the actual Prandtl-Meyer function
  % f_act = f_cal + phi;
  % 
  % % Computing the M_act using Newton-Raphson method
  % tol = 1e-4;
  % error = 1;
  % RF = 0.1;
  % dM = 1e-4;
  % M_act = 1;  % guess value
  % while error > tol
  %   M_act = M_act - RF*(F(f_act, M_act)/F_prime(f_act, M_act, dM));
  %  error = abs(F(f_act, M_act));
  % end
  % 
  % % Computing the actual values of primitive variables on the boundary
  % p_act = p_cal*((1 + ((gamma-1)/2)*M_cal^2)/(1 + ((gamma-1)/2)*M_act^2))^(gamma/(gamma-1));
  % T_act = T_cal*((1 + ((gamma-1)/2)*M_cal^2)/(1 + ((gamma-1)/2)*M_act^2));
  % rho_act = p_act/(R*T_act);
  % 
  % % Updating the boundary node with actual boundary values
  % p(i+1,n_y) = p_act;
  % T(i+1,n_y) = T_act;
  % rho(i+1,n_y) = rho_act;
  % u(i+1,n_y) = u_cal;
  % v(i+1,n_y) = - u_cal*tand(theta);
  % M(i+1,n_y) = M_act;
  % a(i+1,n_y) = (sqrt(gamma*R*T(i+1,n_y)));


  end




















