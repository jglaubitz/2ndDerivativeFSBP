%% script_waveEq
%
% Description: 
%  Script to numerically solve the 1D wave equation equation 
%  Periodic initial and boundary conditions 
%  The FSBP-SAT method is used on a single-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: Feb 14, 2024


%% Setting up the script 
clc, clear 
 

%% Parameters of the problem 
x_L = -1; x_R = 1; % domain boundaries 
c = 1; % speed coefficient 
T = 1; % end time 
% f component 
omega_f = pi; 
% For first test 
f1 = @(x) sin( omega_f*x ); 
f1x = @(x) omega_f*cos( omega_f*x ); 
% For second test 
f2 = @(x) exp( 100*sin( omega_f*x ) ); 
f2x = @(x) 100*omega_f*cos( omega_f*x ).*f2(x); 
% For third test 
f_aux = @(x) x.^2; 
fx_aux = @(x) 2*x;
period = x_R - x_L;
f3 = @(x) f_aux( mod(x + x_L, period) + x_L );
f3x = @(x) fx_aux( mod(x + x_L, period) + x_L );
% g component
omega_g = 2*pi;
g = @(x) cos( omega_g*x ).^2; 
gx = @(x) -2*omega_g*sin( omega_g*x ).*cos( omega_g*x );
% Select components 
f = f3; fx = f3x; % select f component
u_ref = @(x,t) f(x+c*t) + g(x-c*t); % reference solution 
ut_ref = @(x,t) c*fx(x+c*t) - c*gx(x-c*t); % temporal derivative of reference solution
u_init = @(x) u_ref(x,0); % initial data for u 
ut_init = @(x) ut_ref(x,0); % initial data for u_t


%% Shared parameters for the SBP-SAT method 
d = 20;
order = 4; % order of FD-SBP operator  
x_eval = linspace(x_L,x_R,1000); % evaluation points for reference solution
uu_ref = u_ref(x_eval,T); % nodal values of the reference solution 
uut_ref = ut_ref(x_eval,T); % nodal values of the reference solution 


%% Solve the problem using a trigonometric function space on equidistant points 
[D1_trig, D2_trig, x_trig, P_trig, Q_trig ] = compute_FSBP_trig( d ); % construct first- and second-derivative FSBP operators
% Solve the acoustic wave equation 
[ x_trig, u_trig, v_trig, mass_trig, energy_trig ] = solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_trig, D2_trig, x_trig, P_trig );


%% Solve the problem using an FD-SBP operator on the same points  
x_FDSBP = x_trig(1:end-1); % equidistant grid points for FD-SBP operator (exclude right boundary for periodic operator)
N = length(x_FDSBP); % number of equidistant points used by FD-SBP scheme 
dx = (x_R-x_L)/N; % step size between points 
[ D1_FDSBP, D2_FDSBP, P_FDSBP ] = SBP_operators_periodic( N, dx, order ); % get the FD-SBP operators
% Solve the acoustic wave equation 
[ x_FDSBP, u_FDSBP, v_FDSBP, mass_FDSBP, energy_FDSBP ] = solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_FDSBP, D2_FDSBP, x_FDSBP, P_FDSBP ); 


%% Plots 

% Plot for u
figure(1) 
p = plot( x_eval, uu_ref,'k:', x_FDSBP, u_FDSBP,'b--', x_trig, u_trig,'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','FD-SBP','trig. FSBP');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')


% Plot for u_t
figure(2) 
p = plot( x_eval(:), uut_ref(:),'k:', x_FDSBP(:), v_FDSBP(:),'b--', x_trig(:), v_trig(:),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$x$','Interpreter','latex') 
ylabel('$u_t$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','FD-SBP','trig. FSBP');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')


% Plot the mass 
t_aux = mass_trig(:,1); % list of times 
aux = 0;
mass_ref = aux*t_aux.^0; % mass of exact solution over time 
figure(3) 
p = plot( t_aux, mass_ref,'k:', mass_FDSBP(:,1), mass_FDSBP(:,2),'b--', mass_trig(:,1), mass_trig(:,2),'r-.'); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','FD-SBP','trig. FSBP');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')


% Plot the energy 
t_aux = energy_trig(:,1); % list of times 
aux = energy_trig(1,2);
energy_ref = aux*t_aux.^0; % mass of exact solution over time 
figure(4) 
p = plot( t_aux, energy_ref,'k:', energy_FDSBP(:,1), energy_FDSBP(:,2),'b--', energy_trig(:,1), energy_trig(:,2),'r-.' ); 
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('$t$','Interpreter','latex') 
ylabel('$\int u_t^2 + u_x^2 \mathrm{d}x$','Interpreter','latex')
grid on 
lgnd = legend(p, 'ref','FD-SBP','trig. FSBP');
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')