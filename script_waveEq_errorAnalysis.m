%% script_waveEq_errorAnalysis
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
f2 = @(x) exp( 10*sin( omega_f*x ) ); 
f2x = @(x) 10*omega_f*cos( omega_f*x ).*f2(x); 
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
f = f1; fx = f1x; % select f component
u_ref = @(x,t) f(x+c*t) + g(x-c*t); % reference solution 
ut_ref = @(x,t) c*fx(x+c*t) - c*gx(x-c*t); % temporal derivative of reference solution
u_init = @(x) u_ref(x,0); % initial data for u 
ut_init = @(x) ut_ref(x,0); % initial data for u_t


%% Loop over different numbers of points  
d_values = [4,8,16,32,64]; 
N_values = 2*d_values + 2; 
error_trig = []; 
error_FDSBP_order2 = [];
error_FDSBP_order4 = [];
error_FDSBP_order6 = [];
runTime_trig = []; 
runTime_FDSBP_order2 = [];
runTime_FDSBP_order4 = [];
runTime_FDSBP_order6 = [];

for i = 1:length(d_values)

    d = d_values(i)

    %% Construct SBP operators 
    [D1_trig, D2_trig, x_trig, P_trig, Q_trig ] = compute_FSBP_trig( d ); % construct first- and second-derivative FSBP operators
    x_FDSBP = x_trig(1:end-1); % equidistant grid points for FD-SBP method (exclude right boundary for periodic operator)
    N = length(x_FDSBP); % number of equidistant points 
    dx = (x_R-x_L)/N; % spatial step size between points 
    [ D1_FDSBP_order2, D2_FDSBP_order2, P_FDSBP_order2 ] = SBP_operators_periodic( N, dx, 2 ); % get the FD-SBP operators with order 2
    [ D1_FDSBP_order4, D2_FDSBP_order4, P_FDSBP_order4 ] = SBP_operators_periodic( N, dx, 4 ); % get the FD-SBP operators with order 4
    [ D1_FDSBP_order6, D2_FDSBP_order6, P_FDSBP_order6 ] = SBP_operators_periodic( N, dx, 6 ); % get the FD-SBP operators with order 6
    uu_ref_trig = u_ref(x_trig,T); % nodal values of the reference solution
    uu_ref_FDSBP = u_ref(x_FDSBP,T); % nodal values of the reference solution

    % Solve problem with trig FSBP operator and store error 
    [ x, u_trig, v_trig, mass_trig, energy_trig, runTime ] = ... 
        solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_trig, D2_trig, x_trig, P_trig );
    error = sqrt( dot( (u_trig-uu_ref_trig), P_trig*(u_trig-uu_ref_trig) ) ); % absolute error 
    error = error/sqrt( dot( uu_ref_trig, P_trig*uu_ref_trig ) ); % relative error 
    error_trig = [ error_trig; error ]; % save relative error 
    runTime_trig = [ runTime_trig; runTime ]; % save run time 

    % Solve problem with FD-SBP operator with order 2 and store error 
    [ x, u_FDSBP_order2, v_FDSBP_order2, mass_FDSBP_order2, energy_FDSBP_order2, runTime ] = ... 
        solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_FDSBP_order2, D2_FDSBP_order2, x_FDSBP, P_FDSBP_order2 );
    error = sqrt( dot( (u_FDSBP_order2-uu_ref_FDSBP), P_FDSBP_order2*(u_FDSBP_order2-uu_ref_FDSBP) ) ); % absolute error 
    error = error/sqrt( dot( uu_ref_FDSBP, P_FDSBP_order2*uu_ref_FDSBP ) ); % relative error 
    error_FDSBP_order2 = [ error_FDSBP_order2; error ]; % save relative error 
    runTime_FDSBP_order2 = [ runTime_FDSBP_order2; runTime ]; % save run time 

    % Solve problem with FD-SBP operator with order 4 and store error 
    [ x, u_FDSBP_order4, v_FDSBP_order4, mass_FDSBP_order4, energy_FDSBP_order4, runTime ] = ... 
        solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_FDSBP_order4, D2_FDSBP_order4, x_FDSBP, P_FDSBP_order4 );
    error = sqrt( dot( (u_FDSBP_order4-uu_ref_FDSBP), P_FDSBP_order4*(u_FDSBP_order4-uu_ref_FDSBP) ) ); % absolute error 
    error = error/sqrt( dot( uu_ref_FDSBP, P_FDSBP_order4*uu_ref_FDSBP ) ); % relative error 
    error_FDSBP_order4 = [ error_FDSBP_order4; error ]; % save relative error 
    runTime_FDSBP_order4 = [ runTime_FDSBP_order4; runTime ]; % save run time 

    % Solve problem with FD-SBP operator with order 6 and store error 
    [ x, u_FDSBP_order6, v_FDSBP_order6, mass_FDSBP_order6, energy_FDSBP_order6, runTime ] = ... 
        solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1_FDSBP_order6, D2_FDSBP_order6, x_FDSBP, P_FDSBP_order6 );
    error = sqrt( dot( (u_FDSBP_order6-uu_ref_FDSBP), P_FDSBP_order6*(u_FDSBP_order6-uu_ref_FDSBP) ) ); % absolute error 
    error = error/sqrt( dot( uu_ref_FDSBP, P_FDSBP_order6*uu_ref_FDSBP ) ); % relative error 
    error_FDSBP_order6 = [ error_FDSBP_order6; error ]; % save relative error 
    runTime_FDSBP_order6 = [ runTime_FDSBP_order6; runTime ]; % save run time 

end 

%% Plots 

% Plot for relative errors vs number of points N 
figure(1) 
p = plot( N_values, error_FDSBP_order2, 'k--', ...
    N_values, error_FDSBP_order4, 'b:', ... 
    N_values, error_FDSBP_order6, 'g-.', ... 
    N_values, error_trig, 'r-' ...
    ); 
set(p(3), 'color', [0 0.6 0])
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('number of points $N$','Interpreter','latex') 
ylabel('relative $\| \cdot \|_P$-error','Interpreter','latex')
grid on 
lgnd = legend( p, 'FD-SBP, order 2', 'FD-SBP, order 4', 'FD-SBP, order 6', 'trig. FSBP' );
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')

% Normalize the run times
runTime_max = max( [ runTime_FDSBP_order2; runTime_FDSBP_order4; runTime_FDSBP_order6; runTime_trig ] );

% Plot for relative errors vs run time  
figure(2) 
p = plot( runTime_FDSBP_order2/runTime_max, error_FDSBP_order2, 'k--', ...
    runTime_FDSBP_order4/runTime_max, error_FDSBP_order4, 'b:', ... 
    runTime_FDSBP_order6/runTime_max, error_FDSBP_order6, 'g-.', ... 
    runTime_trig, error_trig, 'r-' ... 
    ); 
set(p(3), 'color', [0 0.6 0])
set(p, 'LineWidth',3)
set(gca, 'FontSize', 24)  % Increasing ticks fontsize
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%xlim([x(1),x(end)]) 
%ylim([-1.75,1.75]) 
xlabel('normalized run time','Interpreter','latex') 
ylabel('relative $\| \cdot \|_P$-error','Interpreter','latex')
grid on 
lgnd = legend( p, 'FD-SBP, order 2', 'FD-SBP, order 4', 'FD-SBP, order 6', 'trig. FSBP' );
set(lgnd, 'Interpreter','latex', 'FontSize',24, 'color','none', 'Location','best')