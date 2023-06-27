%% script_viscBurgers_periodic 
%
% Description: 
%  Script to numerically solve the viscous Burgers equation 
%  Periodic initial and boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023

%% Setting up the script 
clc, clear 
 
%% Parameters of the problem 
x_L = 0; x_R = 1; % domain boundaries 
epsilon = 10^(-2); % diffusion parameter 
T = .1; % end time  
u_init = @(x) 2*x + (x >= .5).*( - 2 ); % initial data

%% Shared parameters for the SBP-SAT method 
I = 20; % number of blocks  
I_ref = 200;
d = 2; % degree of the function spaces 
alpha = 1; % parameter of the exponential function space 


%% Solve the problem using a polynomial function space on Lobatto points 
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( d ); % construct first- and second-derivative SBP operators
% Solve problem 
[ x_poly, u_poly, mass_poly, energy_poly ] = solve_viscBurgers_periodic( epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using a exponometric function space on equidistant points 
[D1, D2, x_ref, P, Q] = compute_FSBP_exp( d, alpha, 'LS' ); % construct first- and second-derivative FSBP operators
% Solve problem 
[ x_exp, u_exp, mass_exp, energy_exp ] = solve_viscBurgers_periodic( epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Reference solution using a polynomial function space on Lobatto points with a fine grid
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( d ); % construct first- and second-derivative SBP operators
% Solve problem 
[ x_eval, uu_ref, mass_ref, energy_ref ] = solve_viscBurgers_periodic( epsilon, x_L, x_R, T, u_init, I_ref, D1, D2, x_ref, P );


%% Plot the solutions 
figure(1) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_exp(:), u_exp(:),'r-.' ); 
set(p, 'LineWidth',4)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northeast')

figure(2) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_exp(:), u_exp(:),'r-.' ); 
set(p, 'LineWidth',4)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([.4,.6])  
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northeast')
