%% script_linear_AdvDiff_2D_kernel
%
% Description: 
%  Script to numerically solve the two-dimensional linear advection-diffusion equation 
%  Periodic initial and boundary conditions 
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023 


%% Setting up the script 
%clc, clear 
 

%% Parameters of the problem 
x_L = 0; x_R = 1; % domain boundaries 
a = 1; b = 1; % addvection speeds in x- and y-direction
epsilon = 10^(-4); %10^(-2); % diffusion parameter 
shape_param = 20; % shape parameter of the kernel space 
T = .25; % end time 
u_init = @(x,y) exp(-200*( (x-0.25).^2 + (y-0.25).^2 ) ); % initial condition


%% Shared parameters for the SBP-SAT method  
I = 20; % number of blocks to get numerical solutions 
I_ref = 100; % number of blocks to get reference solution 
x_eval = linspace(x_L,x_R,100); % evaluation points in each direction for reference solution
[xx_init, yy_init] = meshgrid(x_eval, x_eval); 
uu_init = u_init(xx_init,yy_init); % nodal values of the initial data 


%% Solve the problem using a polynomial function space on Lobatto points 
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( 2 ); % construct first- and second-derivative SBP operators
% Solve problem 
[ X_poly, Y_poly, U_poly, mass_poly, energy_poly ] = solve_linear_AdvDiff_2D( a, b, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using a kernel function space on equidistant points (with D1 being correctly (F+F')-exact ) + AD
[D1, D2, x_ref, P, Q] = compute_FSBP_RBF( shape_param, 'LS' ); % construct first- and second-derivative FSBP operators 
% Solve problem 
[ X_RBF, Y_RBF, U_RBF, mass_RBF, energy_RBF ] = solve_linear_AdvDiff_2D( a, b, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Solve the problem using a polynomial function space on equidistant points 
N = length(x_ref); % use the same number of equidistant points as for the exponential function space  
[D1, D2, x_ref, P, Q] = compute_FSBP_poly_equid( 2, N ); % construct first- and second-derivative FSBP operators
% Solve problem 
[ X_poly_equid, Y_poly_equid, U_poly_equid, mass_poly_equid, energy_poly_equid ] = solve_linear_AdvDiff_2D( a, b, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P );


%% Get a reference solution by solving the problem using a polynomial function space on Lobatto points with many blocks  
[D1, D2, x_ref, P, Q] = compute_FSBP_poly( 2 ); % construct first- and second-derivative SBP operators
% Solve problem 
[ X_ref, Y_ref, U_ref, mass_ref, energy_ref ] = solve_linear_AdvDiff_2D( a, b, epsilon, x_L, x_R, T, u_init, I_ref, D1, D2, x_ref, P );


%% store points and the numerical solution a way that allows plotting 

% X-, Y-, and U-values of the numerical solutions  
xx_poly = []; yy_poly = []; uu_poly = []; 
xx_poly_equid = []; yy_poly_equid = []; uu_poly_equid = []; 
xx_RBF = []; yy_RBF = []; uu_RBF = [];
uu_RBF2 = [];
parfor i=1:I 
    xx_poly_aux = []; 
    yy_poly_aux = []; 
    uu_poly_aux = []; 
    xx_poly_equid_aux = []; 
    yy_poly_equid_aux = []; 
    uu_poly_equid_aux = []; 
    xx_RBF_aux = []; 
    yy_RBF_aux = []; 
    uu_RBF_aux = []; 
    uu_RBF2_aux = []; 
    for j=1:I 
        xx_poly_aux = [ xx_poly_aux; X_poly{i,j} ]; 
        yy_poly_aux = [ yy_poly_aux; Y_poly{i,j} ]; 
        uu_poly_aux = [ uu_poly_aux; U_poly{i,j} ]; 
        xx_poly_equid_aux = [ xx_poly_equid_aux; X_poly_equid{i,j} ]; 
        yy_poly_equid_aux = [ yy_poly_equid_aux; Y_poly_equid{i,j} ]; 
        uu_poly_equid_aux = [ uu_poly_equid_aux; U_poly_equid{i,j} ]; 
        xx_RBF_aux = [ xx_RBF_aux; X_RBF{i,j} ]; 
        yy_RBF_aux = [ yy_RBF_aux; Y_RBF{i,j} ]; 
        uu_RBF_aux = [ uu_RBF_aux; U_RBF{i,j} ]; 
        uu_RBF2_aux = [ uu_RBF2_aux; U_RBF{i,j}(1:2:end,1:2:end) ]; 
    end 
    xx_poly = [ xx_poly, xx_poly_aux ];
    yy_poly = [ yy_poly, yy_poly_aux ];
    uu_poly = [ uu_poly, uu_poly_aux ];
    xx_poly_equid = [ xx_poly_equid, xx_poly_equid_aux ];
    yy_poly_equid = [ yy_poly_equid, yy_poly_equid_aux ];
    uu_poly_equid = [ uu_poly_equid, uu_poly_equid_aux ];
    xx_RBF = [ xx_RBF, xx_RBF_aux ];
    yy_RBF = [ yy_RBF, yy_RBF_aux ];
    uu_RBF = [ uu_RBF, uu_RBF_aux ];
    uu_RBF2 = [ uu_RBF2, uu_RBF2_aux ];
end

% X-, Y-, and U-values of the reference solution  
xx_ref = []; yy_ref = []; uu_ref = [];
parfor i=1:I_ref 
    xx_ref_aux = []; 
    yy_ref_aux = []; 
    uu_ref_aux = [];
    for j=1:I_ref 
        xx_ref_aux = [ xx_ref_aux; X_ref{i,j} ]; 
        yy_ref_aux = [ yy_ref_aux; Y_ref{i,j} ]; 
        uu_ref_aux = [ uu_ref_aux; U_ref{i,j} ]; 
    end 
    xx_ref = [ xx_ref, xx_ref_aux ];
    yy_ref = [ yy_ref, yy_ref_aux ];
    uu_ref = [ uu_ref, uu_ref_aux ];
end


%% Compute the pointwise errors 
% We need the reference solution at the same points as the polynomial and
% RBF solution at equidistant grid points 
% To this end, we do a constant interpolation; 

% Make xx_ref, yy_ref into column vectors
xx_ref_vec = xx_ref(:);
yy_ref_vec = yy_ref(:);
uu_ref_vec = uu_ref(:);
xx_poly_vec = xx_poly(:);
yy_poly_vec = yy_poly(:);
uu_poly_vec = uu_poly(:);
uu_RBF_vec = uu_RBF2(:);

% Initialize difference arrays
diff_uu_poly = zeros(size(uu_poly));
diff_uu_RBF = zeros(size(uu_RBF2));

% Loop over points in xx_poly and yy_poly
for i = 1:size(uu_poly, 1)
    for j = 1:size(uu_poly, 2)
        % Find the nearest point in xx_ref, yy_ref to the point in xx_poly, yy_poly
        [~, idx] = min( ( xx_ref_vec - xx_poly(i, j) ).^2 + ( yy_ref_vec - yy_poly(i, j) ).^2 );
        
        % Subtract the corresponding uu_poly value from the nearest uu_ref value
        diff_uu_poly(i, j) = abs( uu_ref_vec(idx) - uu_poly(i, j) ); 

        % Subtract the corresponding uu_RBF2 value from the nearest uu_ref value
        diff_uu_RBF(i, j) = abs( uu_ref_vec(idx) - uu_RBF2(i, j) );

        % Get the value of the reference solution 
        aux(i,j) = uu_ref_vec(idx); 

    end
end

% Compute different relative errors 
error_1Norm_poly = norm( diff_uu_poly(:), 1 ) / norm( aux(:), 1 ); 
error_1Norm_RBF = norm( diff_uu_RBF(:), 1 ) / norm( aux(:), 1 ); 
error_2Norm_poly = norm( diff_uu_poly(:), 2 ) / norm( aux(:), 2 ); 
error_2Norm_RBF = norm( diff_uu_RBF(:), 2 ) / norm( aux(:), 2 ); 
error_maxNorm_poly = norm( diff_uu_poly(:), inf ) / norm( aux(:), inf ); 
error_maxNorm_RBF = norm( diff_uu_RBF(:), inf ) / norm( aux(:), inf ); 


%% Plots 

% Plot the initial data
figure(1) 
s = mesh( xx_init, yy_init, uu_init );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex')

% Plot the initial data
figure(2) 
s = mesh( xx_ref, yy_ref, uu_ref );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex')

% Plot the numerical solution of the polynomial method
figure(3) 
s = mesh( xx_poly, yy_poly, uu_poly );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex')

% Plot the numerical solution of the polynomial method
figure(4) 
s = mesh( xx_poly_equid, yy_poly_equid, uu_poly_equid );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 

% Plot the numerical solution of the kernel method
figure(5) 
s = mesh( xx_RBF, yy_RBF, uu_RBF );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 

% Plot the numerical solution of the kernel method at the same points as
% the polynomial one using Gauss-Lobatto points 
figure(6) 
s = mesh( xx_poly, yy_poly, uu_RBF2 );  
s.EdgeColor = 'interp'; 
set(s, 'LineWidth',1.5); 
set(gca, 'FontSize', 24)  % Increasing ticks fontsize 
xlabel('$x$','Interpreter','latex') 
ylabel('$y$','Interpreter','latex') 
zlabel('$u$','Interpreter','latex') 