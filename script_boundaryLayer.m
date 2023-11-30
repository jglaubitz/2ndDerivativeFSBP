%% script_boundaryLayer
%
% Description: 
%  Script to numerically solve a boundary layer problem for the linear advection-diffusion equation 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023


%% Setting up the script 
clc, clear 
 

%% Parameters of the problem 
x_L = 0; x_R = .5; % domain boundaries 
a = 1; % addvection speed
epsilon = 10^(-2); % diffusion parameter 
T = .75; % end time 
u_init = @(x) 2*x; % initial data 
u_ref = @(x) ( exp(x/epsilon) - 1 ) / ( exp(.5/epsilon) - 1 ); % reference solution 


%% Shared parameters for the SBP-SAT method  
I = 20; % number of blocks  
d = 2; % degree of the function spaces 
alpha = -.1; % parameter of the exponential function space 
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution
uu_ref = u_ref(x_eval); % nodal values of the reference solution 


%% Solve the problem using a polynomial function space on Lobatto points 
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( d ); % construct first- and second-derivative SBP operators
% Solve problem 
[ x_poly, u_poly, mass_poly, energy_poly ] = solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using an exponometric function space on equidistant points 
[D1, D2, x_ref, P, Q] = compute_FSBP_exp( d, alpha, 'LS' ); % construct first- and second-derivative FSBP operators
% Solve problem 
[ x_exp, u_exp, mass_exp, energy_exp ] = solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using a polynomial function space on equidistant points 
N = length(x_ref); % use the same number of equidistant points as for the exponential function space  
[D1, D2, x_ref, P, Q] = compute_FSBP_poly_equid( d, N ); % construct first- and second-derivative FSBP operators
% Solve problem 
[ x_poly_equid, u_poly_equid, mass_poly_equid, energy_poly_equid ] = solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Plots 

% Plot the solutions - Exp vs poly on Guass-Lobatto points  
figure(1) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b-', x_exp(:), u_exp(:),'r--' ); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
ylim([-.1,1]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend( p, 'ref', 'poly', 'exp' );
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northwest')

% Plot the solutions - Exp vs poly on five equidistant points 
figure(2) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly_equid(:), u_poly_equid(:),'g-.', x_exp(:), u_exp(:),'r--' ); 
set(p(2), 'color', [0 0.6 0])
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
ylim([-.1,1]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend( p, 'ref', 'poly', 'exp' );
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northwest')

% Plot the solutions - Zoom-in: Exp vs poly on Guass-Lobatto points 
figure(3) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b-', x_exp(:), u_exp(:),'r-.' ); 
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([0.425,0.5]) 
ylim([-.1,1]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northwest')

% Plot the solutions - Zoom-in: Exp vs poly on five equidistant points 
figure(4) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly_equid(:), u_poly_equid(:),'g-.', x_exp(:), u_exp(:),'r-.' ); 
set(p(2), 'color', [0 0.6 0])
set(p, 'LineWidth',3.5)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([0.425,0.5]) 
ylim([-.1,1]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','exp');
set(lgnd, 'Interpreter','latex', 'FontSize',28, 'color','none', 'Location','northwest')