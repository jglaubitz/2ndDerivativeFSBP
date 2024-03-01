%% solve_waveEq
%
% Description: 
%  Function to solve the wave equation \partial_{tt} u = c^2 \partial_{xx} u 
%  with periodic initial and boundary conditions. 
%  We use a single-block FSBP method to avoid the use of SATs. 
%  (My impression is that finding SATs for the wave equation is a non-trivial and at least partially open problem.) 
%  Time integration is performed with a the SSPRK33 method.  
%
% Author: Jan Glaubitz 
% Date: Feb 15, 2024 
% 
% INPUT: 
%  c :                  speed parameter 
%  x_L, x_R :           left and right boundary of the domain 
%  T :                  end time 
%  u_init, ut_init :    initial data for u and u_t
%  D1, D2 :             PERIODIC first- and second-derivative FSBP operators  
%  x_ref :              grid of operators on reference element [-1,1]
%  P :                  diagonal norm matrix 
%  dt :                 time step size
%
% OUTPUT: 
%  x_global :       global grid points 
%  u_num :          numerical solution at grid points 
%  mass, energy :   mass and energy of the numerical solution over time 
%  runTime :        run time of the simulation 

function [ x_global, u, v, mass, energy, runTime ] = solve_waveEq( c, x_L, x_R, T, u_init, ut_init, D1, D2, x_ref, P, dt )

    %% Set-up the method  
    tic; % start measureing the run time 

    % Data points and the FSBP operators on the physical domain [x_L,x_R]
    N = length(x_ref); % number of data points
    x_global = x_L + (x_R-x_L)*(x_ref + 1)/2; % map points from [-1,1] to [x_L,x_R]
    D1 = (2/(x_R-x_L))*D1; % scale the operator D1
    D2 = (2/(x_R-x_L))^2*D2; % scale the operator D2
    P = .5*(x_R-x_L)*P; % scale the operator P
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 
    
    % Time step size
    dx_min = min(x_global(2:end)-x_global(1:end-1)); % minimum distance between any two neighboring grid points 
    % time step size parameter was not provided, so default is used 
    if ~exist('dt','var')
      dt = 10^(-4); % * dx_min; % time-step size
    end

    % initial data, mass, and energy  
    u = u_init(x_global); % initial values for u on grid   
    v = ut_init(x_global); % initial values for v = u_t on grid   
    mass = []; % mass over time
    energy = []; % energy over time 
     
    % Weak enforcement of boundary conditions and coupling 
    SAT = zeros(N,1); % initialize SAT 
    sigma = 1; 

    %% Iterate over time with the SSPRK33 method until time T is reached 
    t = 0; 
    while (t<T)  

        % time stepping 
        if T-t<dt 
            dt = T-t; 
        else
            t = t+dt;
        end

        % Idea: We re-write the second-order wave equation u_tt = c^2 u_xx as the 
        % first order system u_t = v, v_t = c^2 u_xx with auxilary variable v = u_t. 
        % Then we apply SSPRK33 to that system 
        

        %% SSPRK33
        % 1st update step 
        u_aux1 = u + dt*( v ); % update first component 
        v_aux1 = v + dt*( c^2 * D2*u ); % update second component 
        % 2nd update step 
        u_aux2 = (3/4)*u + (1/4)*u_aux1 + (1/4)*dt*( v_aux1 ); % update first component 
        v_aux2 = (3/4)*v + (1/4)*v_aux1 + (1/4)*dt*( c^2 * D2*u_aux1 ); % update second component 
        % 3th update step 
        u = (1/3)*u + (2/3)*u_aux2 + (2/3)*dt*( v_aux2 ); % update first component 
        v = (1/3)*v + (2/3)*v_aux2 + (2/3)*dt*( c^2 * D2*u_aux2 ); % update second component 

        % Compute mass and energy 
        mass_aux = dot( ones(N,1), P*u ); % approximates \int u(t,x) dx 
        energy_aux = dot( v, P*v ) + dot( D1*u, P*D1*u ); % approximates ||u_t(t,*)||_{L2}^2 + ||u_x(t,*)||_{L2}^2
        % save mass and energy
        mass = [mass; t, mass_aux]; % save mass
        energy = [energy; t, energy_aux]; % save energy 
        
    end
    
    runTime = toc; % stop measuring the run time

end