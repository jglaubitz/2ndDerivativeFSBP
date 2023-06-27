%% compute_FSBP_exp
%
% Description: 
%  Function to compute a first- and second-derivative 
%  exponential FSBP operators on equidistant points on [-1,1]
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT:  
%  d :      degree 
%  alpha :  parameter of the exponential function space 
%  quad :   type of quadrature that is used: LS (least-squares) or GG (generalized Gauss)
%
% OUTPUT: 
%  D1 : differentiation matrix 
%  D2 : differentiation matrix
%  x :  grid points 
%  P :  diagonal-norm matrix
%  Q :  matrix for boundary correction         


function [D1, D2, x, P, Q] = compute_FSBP_exp( d, alpha, quad )

    %% (P1) Compute the first-derivative SBP operator 
    K = d+1; % dimension of the exponential function space  
    
    %% (P2) Generate a basis of F and the corresponding derivatives 
    kk = (0:d-1)'; % exponents 
    basis_F = @(x) [ x.^kk; exp(alpha*x) ]; % basis of F 
    basis_dF = @(x) [ 0*x; kk(2:end).*x.^(kk(2:end)-1); alpha*exp(alpha*x) ]; % the corresponding derivatives 
    
    %% (P3) Compute the quadrature points and weights 
    if strcmp( quad, 'LS')

        % basis of (F^2)', which is needed for LS quadrature 
        kk1 = (0:2*d-3)'; % polynomial degrees 
        kk2 = (0:d-1)'; % degrees of mixed terms poly*exp
        basis_dSqF = @(x) [ x.^(kk1); x.^(kk2).*exp(alpha*x) ; exp(2*alpha*x) ];
        % compute the corresponding moments 
        syms x
        m = int( basis_dSqF(x), [-1 1] ); 
        % compute the points and weights of a positve LS quadrature on
        % [-1,1] with equidistant points 
        [ x, w] = compute_LSQF( -1, 1, basis_dSqF, double(m), 'equid' );

    else

        error('Desired qudarture not yet implemented') 

    end

    %% (P4) Anti-symmetric part of Q 
    N = length(x); % number of grid points
    F = zeros(N,K); % Vandermonde-like matrix 
    dF = zeros(N,K); % Vandermonde-like matrix for the derivatives 
    for n=1:N 
    	F(n,:) = basis_F( x(n) )'; 
        dF(n,:) = basis_dF( x(n) )'; 
    end 
    B = zeros(N); B(1,1) = -1; B(end,end) = 1; % boundary matrix  
    P = diag(w); % diagonal norm matrix contains the quadrature weights 
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