%% solve_viscBurgers_periodic
%
% Description: 
%  Function to numerically solve the viscous Burgers' equation. 
%  Periodic initial and boundary conditions.  
%  The FSBP-SAT method is used on a multi-block structure 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023 
% 
% INPUT: 
%  epsilon :        diffusion parameter 
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

function [ x_global, u, mass, energy ] = solve_linear_AdvDiff( epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P )

    % Data points and the FSBP operator on the reference block [-1,1]
    N = length(x_ref); % number of data points
    block_width = (x_R-x_L)/I; % block width 
    D1 = (2/block_width)*D1; D2 = (2/block_width)^2*D2; P = .5*block_width*P; % scale the operators for local blocks 
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 

    % Global grid points 
    x_global = zeros(N,I); 
    for i=1:I 
        x_global(:,i) = x_L + (x_R-x_L)*(i-1)/I + (x_ref+1)*block_width/2;
    end

    % initial data, mass, and energy  
    u = u_init(x_global); % solution values on the global grid 
    mass = []; % mass over time
    energy = []; % energy over time 

    % Time step size
    dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
    c = max( abs( u(:) ) ); % max absolute advection speed 
    dt = 10^(-2) / ( c*I*N + epsilon*(I*N)^2 ); % time-step size


    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
    t = 0; 
    while (t<T)  

        %% time stepping 
        if T-t<dt 
            dt = T-t; 
        else
            t = t+dt
        end


        %% SSPRK(3,3) time integration 
        % 1st update step 
        SAT = compute_SAT( u, epsilon, D1 );
      	for i = 1:I 
            k1(:,i) = u(:,i) + dt*( -D1*u(:,i).^2/3 - u(:,i).*(D1*u(:,i))/3 + epsilon*D2*u(:,i) + P_inv*SAT(:,i) );  
        end
            
        % 2nd update step 
        SAT = compute_SAT( k1, epsilon, D1 );
        for i = 1:I 
            k1(:,i) = (3/4)*u(:,i) + (1/4)*k1(:,i) + (1/4)*dt*( -D1*k1(:,i).^2/3 - k1(:,i).*(D1*k1(:,i))/3 + epsilon*D2*k1(:,i) + P_inv*SAT(:,i) );  
        end    
        
        % 3th update step 
        SAT = compute_SAT( k1, epsilon, D1 );  
        for i = 1:I
            u(:,i) = (1/3)*u(:,i) + (2/3)*k1(:,i) + (2/3)*dt*( -D1*k1(:,i).^2/3 - k1(:,i).*(D1*k1(:,i))/3 + epsilon*D2*k1(:,i) + P_inv*SAT(:,i) );
        end 
        
        
        %% Compute mass and energy 
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


%% Function to compute BCs and SATs 
function [ SAT ] = compute_SAT( u, epsilon, D1 ) 
    
    % Free parameters of the SATs  
    %sigmaR_1 = a/2; % this choice corresponds to a central flux for advection 
    sigmaR_1 = 0; % this choice corresponds to a full-upwing flux for advection 
    sigmaR_2 = -epsilon/2; % central flux for diffusion 
    
    % All the other parameters follow from these 
    sigmaL_1 = 0; 
    sigmaL_2 = epsilon + sigmaR_2; 
    sigmaL_3 = -sigmaR_2;
    sigmaR_3 = -epsilon - sigmaR_2; 

    % Weak enforcement of boundary conditions and coupling 
    [ N, I ] = size(u);
    SAT_L = zeros(N,1); % initialize left SAT 
    SAT_R = zeros(N,1); % initialize right SAT 
    SAT = zeros(N,I); % initialize overall SAT
    e_L = zeros(N,1); e_L(1) = 1; % auxilary vector used in computing SAT_L 
    e_R = zeros(N,1); e_R(N) = 1; % auxilary vector used in computing SAT_R
    
    for i=1:I 
        
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
        SAT_L = -(2/3)*u_C(1)*e_L*( u_C(1) - u_L(N) ) + sigmaL_2*e_L*( Du_C(1) - Du_L(N) ) + sigmaL_3*(D1.')*e_L*( u_C(1) - u_L(N) ); % left SAT 
        SAT_R = sigmaR_2*e_R*( Du_C(N) - Du_R(1) ) + sigmaR_3*(D1.')*e_R*( u_C(N) - u_R(1) ); % right SAT 
        SAT(:,i) = SAT_L + SAT_R; 
        
    end
    
end