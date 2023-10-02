% Main program file to solve linear convection equation for different
% values of dt

clear all
close all
clc


dt = [1e-1, 1e-2, 1e-3, 1e-4];   % time step

% calling function for solving the linear convection equation and plotting
% the velocity profile for each value of dt

solve_linear_wave_eq_and_plot(dt);

    




