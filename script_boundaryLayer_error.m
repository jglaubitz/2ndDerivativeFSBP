%% script_boundaryLayer_error
%
% Description: 
%  Script to numerically solve a boundary layer problem for the linear advection-diffusion equation 
%  and compare the errors of different methods 
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
d = 2; % degree of the function space 
alpha = -.2; % parameter of the exponential function space 

%% Compute the FSBP operators 

% polynomial function space on Lobatto points 
[D1_poly, D2_poly, x_ref_poly, P_poly, Q] = compute_FSBP_poly( d ); 

% exponometric function space on equidistant points 
[D1_exp, D2_exp, x_ref_exp, P_exp, Q] = compute_FSBP_exp( d, alpha, 'LS' ); 

% polynomial function space on equidistant points 
N = length(x_ref_exp); % use the same number of equidistant points as for the exponential function space  
[D1_poly_equid, D2_poly_equid, x_ref_poly_equid, P_poly_equid, Q] = compute_FSBP_poly_equid( d, N ); % construct first- and second-derivative FSBP operators


%% Prepare loop 
II = [10,20,40,80];
error_L2_poly = zeros(length(II),1); 
error_L2_poly_equid = zeros(length(II),1); 
error_L2_exp = zeros(length(II),1); 
error_max_poly = zeros(length(II),1); 
error_max_poly_equid = zeros(length(II),1); 
error_max_exp = zeros(length(II),1);  

parfor i=1:length(II) 
    
    I = ceil(II(i)) % number of blocks 

    % Solve problem using different FSBP operators
    [ x_poly, u_poly, mass_poly, energy_poly ] = ... 
        solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1_poly, D2_poly, x_ref_poly, P_poly );
    [ x_poly_equid, u_poly_equid, mass_poly_equid, energy_poly_equid ] = ... 
        solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1_poly_equid, D2_poly_equid, x_ref_poly_equid, P_poly_equid );
    [ x_exp, u_exp, mass_exp, energy_exp ] = ... 
        solve_boundaryLayer( a, epsilon, x_L, x_R, T, u_init, I, D1_exp, D2_exp, x_ref_exp, P_exp );

    % compute errors for poly 
    uu_ref = u_ref(x_poly); % values of the steady state solution
    error_L2_poly(i) = norm( uu_ref(:) - u_poly(:), 2 ) / norm( uu_ref(:), 2 ); 
    error_max_poly(i) = norm( uu_ref(:) - u_poly(:), inf ) / norm( uu_ref(:), inf ); 

    % compute errors for poly 
    uu_ref = u_ref(x_poly_equid); % values of the steady state solution
    error_L2_poly_equid(i) = norm( uu_ref(:) - u_poly_equid(:), 2 ) / norm( uu_ref(:), 2 ); 
    error_max_poly_equid(i) = norm( uu_ref(:) - u_poly_equid(:), inf ) / norm( uu_ref(:), inf ); 

    % compute errors for exp 
    uu_ref = u_ref(x_exp); % values of the steady state solution
    error_L2_exp(i) = norm( uu_ref(:) - u_exp(:), 2 ) / norm( uu_ref(:), 2 ); 
    error_max_exp(i) = norm( uu_ref(:) - u_exp(:), inf ) / norm( uu_ref(:), inf ); 

end

%% Plot the errors 

% L2 erros vs I 
figure(1) 
p = plot( II,error_L2_poly,'b^--', II,error_L2_exp,'ro-' ); 
set(p, 'LineWidth',3, 'markersize',14)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xticks(II) 
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_2 / \| u_{\mathrm{ref}} \|_2$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 40*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','Location','best'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 

% Max erros vs I 
figure(2) 
p = plot( II,error_max_poly,'b^--', II,error_max_exp,'ro-' ); 
set(p, 'LineWidth',3, 'markersize',14)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xticks(II) 
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_\infty / \| u_{\mathrm{ref}} \|_\infty$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 40*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','Location','best');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 

% L2 erros vs I 
figure(3) 
p = plot( II,error_L2_poly_equid,'gs-.', II,error_L2_exp,'ro-' ); 
set(p(1), 'color', [0 0.6 0])
set(p, 'LineWidth',3, 'markersize',14)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xticks(II) 
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_2 / \| u_{\mathrm{ref}} \|_2$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log') 
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 40*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','Location','best'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 

% Max erros vs I 
figure(4) 
p = plot( II,error_max_poly_equid,'gs-.', II,error_max_exp,'ro-' );
set(p(1), 'color', [0 0.6 0])
set(p, 'LineWidth',3, 'markersize',14)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
xticks(II) 
xlabel('$I$','Interpreter','latex') 
ylabel('$\| u_{\mathrm{num}} - u_{\mathrm{ref}} \|_\infty / \| u_{\mathrm{ref}} \|_\infty$','Interpreter','latex')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% diagonal line to display rate of convergence 
rate = 2; x = [II(1),II(end)]; y = 40*x.^(-rate); 
rate_line = line(x,y);
rate_line.Color = 'k'; rate_line.LineStyle = ':'; rate_line.LineWidth = 3; 
lgnd = legend('poly','exp','$2$nd order','Location','best'); 
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none')
grid on 