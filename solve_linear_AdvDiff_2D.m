%% solve_linear_AdvDiff_2d
%
% Description: 
%  Function to solve the two-dimensional linear advection-diffusion equation with 
%  periodic initial and boundary conditions using a multi-block FSBP-SAT method. 
%  Time integration with a 3th order TVD/SSP-Runge-Kutta method 
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT: 
%  a, b, epsilon :  advection and diffusion parameters 
%  x_L, x_R :       left and right boundary of the domain 
%  T :              end time 
%  u_init :         initial condition 
%  I :              number of blocks 
%  D1, D2 :         first- and second-derivative FSBP operators  
%  x_ref :          grid of operators on reference element [-1,1]
%  P :              diagonal norm matrix 
%
% OUTPUT: 
%  Cell_X :         Cell structure containing the X-coordinates of the numerical solution 
%  Cell_Y :         Cell structure containing the Y-coordinates of the numerical solution
%  Cell_U :         Cell structure containing the values of the numerical solution
%  mass, energy :   mass and energy of the numerical solution over time    

function [ Cell_X, Cell_Y, Cell_U, mass, energy ] = solve_linear_AdvDiff_2D( a, b, epsilon, x_L, x_R, T, u_init, I, D1, D2, x_ref, P )


    %% Set-up the method  

    % Data points and the FSBP operator on the reference block [-1,1]
    N = length(x_ref); % number of data points
    block_width = (x_R-x_L)/I; % block width 
    D1 = (2/block_width)*D1; % scale the 1st-derivative operator for local blocks 
    D2 = (2/block_width)^2*D2; % scale the 2nd-derivative operator for local blocks
    P = .5*block_width*P; % scale the norm operator for local blocks 
    P_inv = sparse(inv(P)); % precompute inverse diagonal-norm matrix 
    
    % Two-dimensional operators 
    D1x = kron(D1,eye(N)); % 1st-derivative operator in x-direction 
    D1y = kron(eye(N),D1); % 1st-derivative operator in y-direction
    D2x = kron(D2,eye(N)); % 2nd-derivative operator in x-direction 
    D2y = kron(eye(N),D2); % 2nd-derivative operator in y-direction
    P_2d = sparse(kron(P,P)); % two-dimensional norm operator 
    P_2d_inv = inv(P_2d); % inverse two-dimensional norm operator 

    % Time step size
    dx_min = min(x_ref(2:end)-x_ref(1:end-1)); % minimum distance between any two neighboring grid points 
    dt = 10^(-2) / ( max(abs(a),abs(b))*I/dx_min + epsilon*(I/dx_min)^2 ); % time-step size

    % Global grid points in each direction 
    x_1d = zeros(N,I); 
    for i=1:I 
        x_1d(:,i) = x_L + (x_R-x_L)*(i-1)/I + (x_ref+1)*block_width/2;
    end
    
    % Set-up a cell structure that contains the two-dimensional grid points 
    Cell_X = cell(I,I); % x-coordinates in each block 
    Cell_Y = cell(I,I); % y-coordinates in each block 
    Cell_U = cell(I,I); % nodal values of the numerical solution in each block 
    Cell_K1 = cell(I,I); % cell axilary matrices to update the numerical solution in each block 
    for i=1:I % loop over blocks in x-direction 
        for j=1:I % loop over blocks in y-direction 
            [xx, yy] = meshgrid(x_1d(:,i), x_1d(:,j)); % generate 2d points 
            Cell_X{i,j} = xx; % store x-coordinates 
            Cell_Y{i,j} = yy; % store y-coordinates 
            Cell_U{i,j} = u_init(xx,yy); % store store nodal values of initial data 
        end
    end

    % initiate mass and energy  
    mass = []; % mass over time
    energy = []; % energy over time 
     
    
    %% Iterate over time with a 3th-order Runge-Kutta until time T is reached 
    t = 0; 
    counter = 0; 
    while (t<T)  

        % time stepping 
        if T-t<dt 
            dt = T-t; 
        else
            t = t+dt
        end

        % SSPRK(3,3) time integration 
        % 1st update step 
        SAT = compute_SAT( Cell_U, a, b, epsilon, D1, I );
      	parfor i = 1:I 
            for j = 1:I 
                % vectorize solution and SAT matrix  
                u = reshape( Cell_U{i,j}, N^2, 1 ); 
                sat = reshape( SAT{i,j}, N^2, 1 );
                % update solution 
                k1 = u + dt*( ... 
                    -( a*D1x + b*D1y )*u + ... 
                    epsilon*( D2x + D2y )*u + ... 
                    P_2d_inv*sat ); 
                % Transform into matrix again 
                Cell_K1{i,j} = reshape( k1, N, N ); 
            end
        end
            
        % 2nd update step 
        SAT = compute_SAT( Cell_K1, a, b, epsilon, D1, I );
        parfor i = 1:I 
            for j = 1:I
                % vectorize solution and SAT matrices  
                u = reshape( Cell_U{i,j}, N^2, 1 ); 
                k1 = reshape( Cell_K1{i,j}, N^2, 1 ); 
                sat = reshape( SAT{i,j}, N^2, 1 );
                % update solution 
                k1 = (3/4)*u + (1/4)*k1 + (1/4)*dt*( ... 
                    -( a*D1x + b*D1y )*k1 + ... 
                    epsilon*( D2x + D2y )*k1 + ... 
                    P_2d_inv*sat ); 
                % Transform into matrix again 
                Cell_K1{i,j} = reshape( k1, N, N ); 
            end    
        end

        % 3th update step 
        SAT = compute_SAT( Cell_K1, a, b, epsilon, D1, I ); 
        parfor i = 1:I 
            for j = 1:I 
                % vectorize solution and SAT matrices  
                u = reshape( Cell_U{i,j}, N^2, 1 ); 
                k1 = reshape( Cell_K1{i,j}, N^2, 1 ); 
                sat = reshape( SAT{i,j}, N^2, 1 );
                % update solution 
                u = (1/3)*u + (2/3)*k1 + (2/3)*dt*( ... 
                    -( a*D1x + b*D1y )*k1 + ... 
                    epsilon*( D2x + D2y )*k1 + ... 
                    P_2d_inv*sat ); 
                % Transform into matrix again 
                Cell_U{i,j} = reshape( u, N, N ); 
            end
        end 
        
        
        %% Compute mass and energy 
        counter = counter + 1; 
        if mod(counter,10) == 1
            mass_aux = 0; energy_aux = 0; 
           parfor i=1:I 
               for j=1:I
                   mass_aux = mass_aux + dot( diag(P_2d), reshape( Cell_U{i,j}, N^2, 1 ) ); % compute mass 
                   energy_aux = energy_aux + dot( diag(P_2d), reshape( Cell_U{i,j}.^2, N^2, 1 ) ); % compute energy
               end
           end 
            % save mass and energy
            mass = [mass; t, mass_aux]; % save mass
            energy = [energy; t, energy_aux]; % save energy
        end 

    end
    
end


%% Function to compute BCs and SATs 
function [ SAT ] = compute_SAT( U, a, b, epsilon, D1, I ) 
    
    % Free parameters of the SATs  
    sigmaR_1 = .5*a; % this choice corresponds to a central flux for advection 
    %sigmaR_1 = 0; % this choice corresponds to a full-upwing flux for advection 
    sigmaR_2 = -epsilon/2; % central flux for diffusion 
    
    sigmaT_1 = .5*b; % this choice corresponds to a central flux for advection 
    %sigmaT_1 = 0; % this choice corresponds to a full-upwing flux for advection 
    sigmaT_2 = -epsilon/2; % central flux for diffusion 
    
    % All the other parameters follow from these 
    sigmaL_1 = sigmaR_1 - a; 
    sigmaL_2 = epsilon + sigmaR_2; 
    sigmaL_3 = -sigmaR_2;
    sigmaR_3 = -epsilon - sigmaR_2; 

    sigmaB_1 = sigmaT_1 - b; 
    sigmaB_2 = epsilon + sigmaT_2; 
    sigmaB_3 = -sigmaT_2;
    sigmaT_3 = -epsilon - sigmaT_2; 

    % Weak enforcement of boundary conditions and coupling  
    SAT = cell(I,I); % initialize a cell structure with a SAT for each block 
    N = size(D1,1); % number of grid points in each direction 
    D1x = kron(D1,eye(N)); % 1st-derivative operator in x-direction 
    D1y = kron(eye(N),D1); % 1st-derivative operator in y-direction

    parfor i=1:I % loop over blocks in x-direction 
        for j=1:I % loop over blocks in y-direction 
        
            %% initialize SATs in different direction 
            SAT_L = zeros(N,N); % initialize left SAT 
            SAT_R = zeros(N,N); % initialize right SAT 
            SAT_B = zeros(N,N); % initialize bottom SAT 
            SAT_T = zeros(N,N); % initialize top SAT
 
            %% get the solution in the block at hand and the neighboring blocks 
            %% We assume periodic BCs 

            % solution in the block at hand 
            U_C = U{i,j}; 
            
            % solution in the block to the left-hand side 
            if i==1 % if we are considering the first block in a row  
                U_L = U{I,j}; % we take the last block in this row 
            else % otherwise 
                U_L = U{i-1,j}; % we take the block to the left-hand side
            end

            % solution in the block to the right-hand side 
            if i==I % if we are considering the last block in a row  
                U_R = U{1,j}; % we take the first block in this row
            else % otherwise 
                U_R = U{i+1,j}; % we take the block to the right-hand side
            end

            % solution in the block below 
            if j==1 % if we are considering the first block in a row  
                U_B = U{i,I}; % we take the last block in this row 
            else % otherwise 
                U_B = U{i,j-1}; % we take the block to the left-hand side
            end

            % solution in the block on top 
            if j==I % if we are considering the last block in a row  
                U_T = U{i,1}; % we take the first block in this row
            else % otherwise 
                U_T = U{i,j+1}; % we take the block to the right-hand side
            end

            % Vectorize the solution matrices 
            u_C = U_C(:); 
            u_L = U_L(:); 
            u_R = U_R(:); 
            u_B = U_B(:); 
            u_T = U_T(:); 

            % compute the derivatives 
            Dx_u_C = D1x*u_C; 
            Dy_u_C = D1y*u_C; 
            Dx_u_L = D1x*u_L; 
            Dy_u_L = D1y*u_L;
            Dx_u_R = D1x*u_R; 
            Dy_u_R = D1y*u_R; 
            Dx_u_B = D1x*u_B; 
            Dy_u_B = D1y*u_B;
            Dx_u_T = D1x*u_T; 
            Dy_u_T = D1y*u_T;
            
            % re-write as matrices 
            Dx_U_C = reshape( Dx_u_C, N, N );  
            Dy_U_C = reshape( Dy_u_C, N, N ); 
            Dx_U_L = reshape( Dx_u_L, N, N );  
            Dy_U_L = reshape( Dy_u_L, N, N ); 
            Dx_U_R = reshape( Dx_u_R, N, N );  
            Dy_U_R = reshape( Dy_u_R, N, N );  
            Dx_U_B = reshape( Dx_u_B, N, N ); 
            Dy_U_B = reshape( Dy_u_B, N, N ); 
            Dx_U_T = reshape( Dx_u_T, N, N );  
            Dy_U_T = reshape( Dy_u_T, N, N ); 
            
            Dx_diff_L = zeros(N,N); 
            Dx_diff_L(:,1) = U_C(:,1) - U_L(:,end); 
            DxT_L = reshape( (D1x.')*Dx_diff_L(:), N, N ); 
            
            Dx_diff_R = zeros(N,N); 
            Dx_diff_R(:,end) = U_C(:,end) - U_R(:,1); 
            DxT_R = reshape( (D1x.')*Dx_diff_R(:), N, N );

            Dx_diff_B = zeros(N,N); 
            Dx_diff_B(1,:) = U_C(1,:) - U_B(end,:); 
            DxT_B = reshape( (D1y.')*Dx_diff_B(:), N, N );

            Dx_diff_T = zeros(N,N); 
            Dx_diff_T(end,:) = U_C(end,:) - U_T(1,:); 
            DxT_T = reshape( (D1y.')*Dx_diff_T(:), N, N );

            %% get the SATs 
            % left SAT 
            SAT_L(:,1) = sigmaL_1*( U_C(:,1) - U_L(:,end) ) + ... 
                sigmaL_2*( Dx_U_C(:,1) - Dx_U_L(:,end) ); 
            %SAT_L = SAT_L + sigmaL_3*DxT_L; 
            % right SAT
            SAT_R(:,end) = sigmaR_1*( U_C(:,end) - U_R(:,1) ) + ... 
                sigmaR_2*( Dx_U_C(:,end) - Dx_U_R(:,1) );  
            %SAT_R = SAT_R + sigmaR_3*DxT_R; 
            % bottom SAT 
            SAT_B(1,:) = sigmaB_1*( U_C(1,:) - U_B(end,:) ) + ... 
                sigmaB_2*( Dy_U_C(1,:) - Dy_U_B(end,:) );  
            %SAT_B = SAT_B + sigmaB_3*DxT_B; 
            % top SAT 
            SAT_T(end,:) = sigmaT_1*( U_C(end,:) - U_T(1,:) ) + ... 
                sigmaT_2*( Dy_U_C(end,:) - Dy_U_T(1,:) ); 
            %SAT_T = SAT_T + sigmaT_3*DxT_T; 
            % sum-up SATs
            SAT{i,j} = SAT_L + SAT_R + SAT_B + SAT_T; 

        end
    end
    
end