%% script_waveEq_specturm
%
% Description: 
%  Script to numerically solve the 1D wave equation equation 
%  Periodic initial and boundary conditions 
%  The FSBP-SAT method is used on a single-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: Feb 16, 2024


%% Setting up the script 
clc, clear 
x_L = -1; x_R = 1;
d = 20; 


%% Trigonometric FSBP operator on equidistant points 
[D1_trig, D2_trig, x_trig, P_trig, Q_trig ] = compute_FSBP_trig( d ); % construct first- and second-derivative FSBP operators
lambda_D1_trig = eig(D1_trig); % eigenvalues of first-derivative operator 
lambda_D2_trig = eig(D2_trig); % eigenvalues of second-derivative operator


%% Traditional FD-SBP operator of order 2 on the same points  
x_FDSBP = x_trig; % equidistant grid points for FD-SBP operator
N = length(x_FDSBP); % number of equidistant points 
dx = (x_R-x_L)/(N-1); % step size 
[ D1_FDSBP_order2, D2_FDSBP_order2, P_FDSBP_order2 ] = SBP_operators_periodic( N, dx, 2 ); % get the FD-SBP operators
lambda_D1_FDSBP_order2 = eig(D1_FDSBP_order2); % eigenvalues of first-derivative operator 
lambda_D2_FDSBP_order2 = eig(D2_FDSBP_order2); % eigenvalues of second-derivative operator


%% Traditional FD-SBP operator of order 4 on the same points  
[ D1_FDSBP_order4, D2_FDSBP_order4, P_FDSBP_order4 ] = SBP_operators_periodic( N, dx, 4 ); % get the FD-SBP operators
lambda_D1_FDSBP_order4 = eig(D1_FDSBP_order4); % eigenvalues of first-derivative operator 
lambda_D2_FDSBP_order4 = eig(D2_FDSBP_order4); % eigenvalues of second-derivative operator


%% Traditional FD-SBP operator of order 6 on the same points  
[ D1_FDSBP_order6, D2_FDSBP_order6, P_FDSBP_order6 ] = SBP_operators_periodic( N, dx, 6 ); % get the FD-SBP operators
lambda_D1_FDSBP_order6 = eig(D1_FDSBP_order6); % eigenvalues of first-derivative operator 
lambda_D2_FDSBP_order6 = eig(D2_FDSBP_order6); % eigenvalues of second-derivative operator


%% Plot the results 

% Spectrum of D2: trig FSBP vs FD-SBP with order 2
figure(1) 
p1 = plot( real(lambda_D2_trig), imag(lambda_D2_trig), 'ro', ... 
    real(lambda_D2_FDSBP_order2), imag(lambda_D2_FDSBP_order2), 'b*' ); 
set(p1, 'markersize',14, 'LineWidth',3); 
set(gca, 'FontSize', 24); % Increasing ticks fontsize 
xlabel('real part','Interpreter','latex'); 
ylabel('imaginary part','Interpreter','latex'); 
grid on 
%lgnd = legend('true signal $x$','noisy data $\mathbf{y}$'); 
lgnd = legend('trig FSBP','FD-SBP, order two');  
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none', 'Location','northwest')
hold off

% Spectrum of D2: trig FSBP vs FD-SBP with order 4
figure(2) 
p1 = plot( real(lambda_D2_trig), imag(lambda_D2_trig), 'ro', ... 
    real(lambda_D2_FDSBP_order4), imag(lambda_D2_FDSBP_order4), 'b*' ); 
set(p1, 'markersize',14, 'LineWidth',3); 
set(gca, 'FontSize', 24); % Increasing ticks fontsize 
xlabel('real part','Interpreter','latex'); 
ylabel('imaginary part','Interpreter','latex'); 
grid on 
%lgnd = legend('true signal $x$','noisy data $\mathbf{y}$'); 
lgnd = legend('trig FSBP','FD-SBP, order four');  
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none', 'Location','northwest')
hold off

% Spectrum of D2: trig FSBP vs FD-SBP with order 6
figure(3) 
p1 = plot( real(lambda_D2_trig), imag(lambda_D2_trig), 'ro', ... 
    real(lambda_D2_FDSBP_order6), imag(lambda_D2_FDSBP_order6), 'b*' ); 
set(p1, 'markersize',14, 'LineWidth',3); 
set(gca, 'FontSize', 24); % Increasing ticks fontsize 
xlabel('real part','Interpreter','latex'); 
ylabel('imaginary part','Interpreter','latex'); 
grid on 
%lgnd = legend('true signal $x$','noisy data $\mathbf{y}$'); 
lgnd = legend('trig FSBP','FD-SBP, order six');  
set(lgnd, 'Interpreter','latex', 'FontSize',26, 'color','none', 'Location','northwest')
hold off