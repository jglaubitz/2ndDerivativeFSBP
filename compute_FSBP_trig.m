%% compute_FSBP_trig
%
% Description: 
%  Function to compute a first- and second-derivative 
%  triginometric FSBP operators on equidistant points on [-1,1]
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT:  
%  d :  degree 
%
% OUTPUT: 
%  D1 : differentiation matrix 
%  D2 : differentiation matrix
%  x :  grid points 
%  P :  diagonal-norm matrix
%  Q :  matrix for boundary correction         

function [D1, D2, x, P, Q ] = compute_FSBP_trig( d )

    % Make sure that d is large enough 
    if d < 1 
        error('degree d too small')
    end

    %% (P1) Compute the first-derivative SBP operator 
    K = 2*d+1; % dimension of the trigonometric function space 
    N = K+1; % number of grid points (need one more for unisolvency) 
    
    %% (P2) Compute the trapezoidal rule  
   	x = linspace(-1, 1, N)'; 
   	w = 2/(N-1)*ones(N,1); 
   	w(1) = 0.5*w(1); w(end) = 0.5*w(end);  
    P = diag(w); % diagonal norm matrix 

    %% (P3) Generate a basis of F and F' 
    kk = (1:d)'; 
    % basis of F and the corresponding derivatives 
    basis_F = @(x) [x.^0; sin( kk*pi*x ); cos( kk*pi*x ) ];
    basis_dF = @(x) [0*x; kk*pi.*cos( kk*pi*x ); -kk*pi.*sin( kk*pi*x ) ];

    %% (P4) Anti-symmetric part of Q 
    F = zeros(N,K); % Vandermonde-like matrix 
    dF = zeros(N,K); % Vandermonde-like matrix for the derivatives 
    for n=1:N 
    	F(n,:) = basis_F( x(n) )'; 
        dF(n,:) = basis_dF( x(n) )'; 
    end 
    B = zeros(N); B(1,1) = -1; B(end,end) = 1; % boundary matrix  
    RHS = P*dF - 0.5*B*F; % right-hand side of the matrix equation for Q_anti 
    
    % Vectorize 
    A = kron( F', eye(N) ); % coefficient matrix for vectorized version 
    rhs = RHS(:); % vectorized right-hand side 
    
    % Commutation matrix C
    I = reshape(1:N*N, [N, N]); % initialize a matrix of indices 
    I = I'; % transpose it
    I = I(:); % vectorize the required indices
    C = speye(N*N); % Initialize an identity matrix
    C = C(I,:); % Re-arrange the rows of the identity matrix
    
    % Get the anti-symmetric least-squares solution to 
    A_ext = [ A; C+speye(N^2) ]; 
    rhs_ext = [ rhs; zeros(N^2,1) ]; 
    q_anti = lsqminnorm(A_ext,rhs_ext); % least-squares solution  
    Q_anti = reshape(q_anti,N,N);

    %% (P5) Compute Q, D1, and D2
    Q = Q_anti + 0.5*B; % matrix Q 
    D1 = diag(1./w)*Q; % first-derivative FSBP operator 
    D2 = D1*D1; % second-derivative FSBP operator 
    
end