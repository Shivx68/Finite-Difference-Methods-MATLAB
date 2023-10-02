% Main script- Simulation of quasi 1D subsonic-supersonic isentropic flow 
% nozzle in non-conservation and conservation forms

clear all
close all
clc

L = 3;
n = 61;
throat = (n+1)/2;
x = linspace(0,3,n);
nt = 1400;
C = 0.5;
tol = 1e-4;
m_analytical = 0.579*ones(n);


##% solving governing equations in the non-conservation form
[rho_nc,T_nc,p_nc,M_nc,m_nc,n_nc] = non_conservation_form(L,n,nt,C,tol);

% solving giverning equations in the consevtion form
%[rho_c,T_c,p_c,M_c,m_c,n_c] = conservation_form(L,n,nt,C,tol);

##% Comparing the normalized mass flow rate between conservation &
##% non-conservation forms
##figure(7)
##plot(x,m_nc,x,m_c,x,m_analytical,'g','LineWidth',2);
##xlabel('x/L');
##ylabel('\rhoVA/\rho_oV_oA*');
##title(sprintf(['Comparison of steady state mass flow variations obtained with'...
##    '\nconservation and non-conservation forms']));
##legend('Non-conservation form','Conservation form','Analytical solution');
##
##% comparing the convergence rate
##figure(8)
##c = categorical({'Non-conservation form','Conservation form'});
##b = bar(c, [n_nc n_c]);
##b.FaceColor = 'flat';
##b.CData(2,:) = [1 1 0];
##ylabel('No. of iterations');
##title(sprintf(['Numer of iterations taken to converge to a steady-state'...
##    '\n for a tolerance of %g'], tol));
##
##% Computing execution times accurately using timeit()
##myfun_nc = @()non_conservation_form(L,n,nt,C,tol);
##time_nc = timeit(myfun_nc,6);
##myfun_c = @()conservation_form(L,n,nt,C,tol);
##time_c = timeit(myfun_c,6);
##
##% displaying the execution times on the output window
##fprintf('Time taken to execute non-conservation form = %0.3gs', time_nc);
##fprintf('\nTime taken to execute conservation form = %0.3gs\n', time_c);
%Calculations for graphic display start

