%% solve_linear_AdvDiff
%
% Description: 
%  Function to solve the linear advection-diffusion equation with periodic 
%  initial and boundary conditions using a single/multi-block FSBP-SAT method. 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023 
% 
% INPUT: 
%  a, epsilon :     advection and diffusion parameter 
%  x_L, x_R :       left and right boundary of the domain 
%  T :              end time 
%  u_init :         initial condition 
%  I :              number of blocks 
%  D1, D2 :         first- and second-derivative FSBP operators  
%  x_ref :          grid of operators on reference element [-1,1]
%  P :              diagonal norm matrix 
%
% OUTPUT: 
%  x_global :       global grid points 
%  u_num :          numerical solution at grid points 
%  mass, energy :   mass and energy of the numerical solution over time    

function [ x_global, u_num, mass, energy ] = solve_linear_AdvDiff( a, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P )

    %% Set-up the method  

    % Data points and the FSBP operator on the reference block [-1,1]
    N = length(x_ref); % number of data points
    block_width = (x_R-x_L)/I; % block width 
    D1 = (2/block_width)*D1; D2 = (2/block_width)^2*D2; P = .5*block_width*P; % scale the operators for local blocks 
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 
    
    % Time step size
    dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
    dt = 10^(-2) / ( a*I/dx_min + epsilon*(I/dx_min)^2 ); % time-step size

    % Global grid points 
    x_global = zeros(N,I); 
    for i=1:I 
        x_global(:,i) = x_L + (x_R-x_L)*(i-1)/I + (x_ref+1)*block_width/2;
    end

    % initial data, mass, and energy  
    u = u_init(x_global); % solution values on the global grid 
    mass = []; % mass over time
    energy = []; % energy over time 


    %% Set-up the SATs for 

    % Free parameters of the SATs  
    %sigmaR_1 = a/2; % this choice corresponds to a central flux for advection 
    sigmaR_1 = 0; % this choice corresponds to a full-upwing flux for advection 
    sigmaR_2 = -epsilon/2; % central flux for diffusion 
    
    % All the other parameters follow from these 
    sigmaL_1 = sigmaR_1 - a; 
    sigmaL_2 = epsilon + sigmaR_2; 
    sigmaL_3 = -sigmaR_2;
    sigmaR_3 = -epsilon - sigmaR_2; 

    % Weak enforcement of boundary conditions and coupling 
    SAT_L = zeros(N,1); % initialize left SAT 
    SAT_R = zeros(N,1); % initialize right SAT 
    SAT = zeros(N,1); % initialize overall SAT
    e_L = zeros(N,1); e_L(1) = 1; % auxilary vector used in computing SAT_L 
    e_R = zeros(N,1); e_R(N) = 1; % auxilary vector used in computing SAT_R 
     
    
    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
    t = 0; 
    while (t<T)  

        % time stepping 
        if T-t<dt 
            dt = T-t; 
        else
            t = t+dt
        end

        % loop over block 
        for i = 1:I

            %% SATs to weakly enforce BC and coupling 
            u_C = u(:,i); % solution in the element at hand 
            
            % solution in the element to the left-hand side 
            if i==1 % if we are considering the first element 
                u_L = u(:,I); % we take the last element, assuming periodic BC
            else % otherwise 
                u_L = u(:,i-1); % we take the element to the left-hand side
            end

            % solution in the element to the right-hand side 
            if i==I % if we are considering the last element 
                u_R = u(:,1); % we take the first element, assuming periodic BC
            else % otherwise 
                u_R = u(:,i+1); % we take the element to the right-hand side
            end

            % compute the derivatives of u_L, u_C, and u_R 
            Du_L = D1*u_L; 
            Du_C = D1*u_C; 
            Du_R = D1*u_R;

            % get the SATs 
            SAT_L = sigmaL_1*e_L*( u_C(1) - u_L(N) ) + sigmaL_2*e_L*( Du_C(1) - Du_L(N) ) + sigmaL_3*(D1.')*e_L*( u_C(1) - u_L(N) ); % left SAT 
            SAT_R = sigmaR_1*e_R*( u_C(N) - u_R(1) ) + sigmaR_2*e_R*( Du_C(N) - Du_R(1) ) + sigmaR_3*(D1.')*e_R*( u_C(N) - u_R(1) ); % right SAT 
            SAT = SAT_L + SAT_R; 


            %% FSBP
            % 1st update step 
            k1 = u(:,i) + dt*( -a*D1*u(:,i) + epsilon*D2*u(:,i) + P_inv*SAT );  
            % 2nd update step 
            k1 = (3/4)*u(:,i) + (1/4)*k1 + (1/4)*dt*( -a*D1*k1 + epsilon*D2*k1 + P_inv*SAT );  
            % 3th update step 
            u_num(:,i) = (1/3)*u(:,i) + (2/3)*k1 + (2/3)*dt*( -a*D1*k1 + epsilon*D2*k1 + P_inv*SAT );

        end

        u = u_num; % update solution 
        
        % Compute mass and energy 
        mass_aux = 0; energy_aux = 0; 
        for i=1:I 
            mass_aux = mass_aux + dot( ones(N,1), P*u(:,i) ); % compute mass 
            energy_aux = energy_aux + dot( u(:,i), P*u(:,i) ); % compute energy
        end 
        % save mass and energy
        mass = [mass; t, mass_aux]; % save mass
        energy = [energy; t, energy_aux]; % save energy 
        
    end
    
end