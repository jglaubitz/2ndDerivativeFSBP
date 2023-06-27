%% compute_FSBP_RBF
%
% Description: 
%  Function to compute a first- and second-derivative 
%  RBF FSBP operators on equidistant points on [-1,1]
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT:  
%  shape_param :    shape parameter of the Gaussian RBF 
%  quad :           type of quadrature that is used: LS (least-squares) or GG (generalized Gauss)
%
% OUTPUT: 
%  D1 : differentiation matrix 
%  D2 : differentiation matrix
%  x :  grid points 
%  P :  diagonal-norm matrix
%  Q :  matrix for boundary correction         

function [D1, D2, x, P, Q] = compute_FSBP_RBF( shape_param, quad )

    %% (P1) Compute the first-derivative SBP operator 
    K = 4; % dimension of the RBF function space  
    
    %% (P2) Generate a basis of G = F + F' and the corresponding derivatives 
    basis_G = @(x) [ x.^0; x; exp( -(x/shape_param).^2 ); x.*exp( -(x/shape_param).^2 ) ]; % basis of G = F + F'
    basis_dG = @(x) [ 0*x; x.^0; -2*x/(shape_param^2)*exp( -(x/shape_param).^2 ) ; ... 
        (1 -2*(x/shape_param).^2)*exp( -(x/shape_param).^2 ) ]; % the corresponding derivatives 

    %% (P3) Compute the quadrature points and weights 
    if strcmp( quad, 'LS')

        % basis of (F^2)', which is needed for LS quadrature 
        basis_dSqG = @(x) [ ...
            x.^0; x; ... 
            x.*exp( -(x/shape_param).^2 ); ... 
            ( 1 - 2*(x/shape_param).^2 ).*exp( -(x/shape_param).^2 ); ...
            2*x.*( 1 - (x/shape_param).^2 ).*exp( -(x/shape_param).^2 ); ...
            x.*exp( -2*(x/shape_param).^2 ); ...
            ( 1 - 4*(x/shape_param).^2 ).*exp( -2*(x/shape_param).^2 ); ...
            2*x.*( 1 - 2*(x/shape_param).^2 ).*exp( -2*(x/shape_param).^2 ); ...
            ];
        % compute the corresponding moments 
        syms x
        m = int( basis_dSqG(x), [-1 1]); 
        % compute the points and weights of a positve LS quadrature on
        % [-1,1] with equidistant points 
        [ x, w] = compute_LSQF( -1, 1, basis_dSqG, double(m), 'equid' );

    else 

        error('Desired qudarture not yet implemented') 

    end

    %% (P4) Anti-symmetric part of Q 
    N = length(x); % number of grid points
    G = zeros(N,K); % Vandermonde-like matrix 
    dG = zeros(N,K); % Vandermonde-like matrix for the derivatives 
    for n=1:N 
    	G(n,:) = basis_G( x(n) )'; 
        dG(n,:) = basis_dG( x(n) )'; 
    end 
    B = zeros(N); B(1,1) = -1; B(end,end) = 1; % boundary matrix  
    P = diag(w); % diagonal norm matrix contains the quadrature weights 
    RHS = P*dG - 0.5*B*G; % right-hand side of the matrix equation for Q_anti 
    
    % Vectorize 
    A = kron( G', eye(N) ); % coefficient matrix for vectorized version 
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