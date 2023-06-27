%% script_linear_AdvDiff
%
% Description: 
%  Script to numerically solve the linear advection-diffusion equation 
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
a = 1; % addvection speed
epsilon = 10^(-5); % diffusion parameter 
T = 1; % end time 
omega1 = 4*pi; omega2 = 40*pi; % frequencies of the initial data 
u_init = @(x) cos(omega1*x) + 0.75*sin(omega2*x); % initial data 
u_ref = @(x,t) exp( -epsilon*(omega1^2)*t ).*cos(omega1*(x-t)) + ... 
    0.75*exp( -epsilon*(omega2^2)*t ).*sin(omega2*(x-t)); % reference solution 


%% Shared parameters for the SBP-SAT method  
I = 1; % number of blocks  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution
uu_ref = u_ref(x_eval,T); % nodal values of the reference solution 


%% Solve the problem using a polynomial function space on Lobatto points 
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( 60 ); % construct first- and second-derivative SBP operators
% Solve problem 
[ x_poly, u_poly, mass_poly, energy_poly ] = solve_linear_AdvDiff( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using a trigonometric function space on equidistant points 
[D1, D2, x_ref, P, Q] = compute_FSBP_trig( 30 ); % construct first- and second-derivative FSBP operators
% Solve problem 
[ x_trig, u_trig, mass_trig, energy_trig ] = solve_linear_AdvDiff( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Plots 

% Plot the solutions 
figure(1) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_trig(:), u_trig(:),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
ylim([-1.75,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Plot the zoomed-in solutions 
figure(2) 
p = plot( x_eval(:), uu_ref(:),'k:', x_poly(:), u_poly(:),'b--', x_trig(:), u_trig(:),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xlim([.4,.6]) 
ylim([-.25,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')


% Plot the mass 
t_aux = mass_trig(:,1); % list of times 
aux = 0;
mass_ref = aux*t_aux.^0; % mass of exact solution over time 
figure(3) 
p = plot( t_aux, mass_ref,'k:', mass_poly(:,1), mass_poly(:,2),'b--', mass_trig(:,1), mass_trig(:,2),'r-.'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Compute the reference solution's energy over time 
t_aux = energy_trig(:,1); % list of times 
alpha = exp( -epsilon*(omega1^2)*t_aux ); % time-dependent coefficient
beta = 0.75*exp( -epsilon*(omega2^2)*t_aux ); % time-dependent coefficient 
syms x; 
int_aux1 = int( cos(omega1*x).^2, x_L,x_R);
int_aux2 = int( cos(omega1*x).*sin(omega2*x), x_L,x_R);
int_aux3 = int( sin(omega2*x).^2, x_L,x_R);
energy_ref = (alpha.^2)*int_aux1 + 2*alpha.*beta*int_aux2 + (beta.^2)*int_aux3;

% Plot the energy 
figure(4) 
p = plot( t_aux, energy_ref,'k:', energy_poly(:,1), energy_poly(:,2),'b--', energy_trig(:,1), energy_trig(:,2),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u^2 \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','poly','trig');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')