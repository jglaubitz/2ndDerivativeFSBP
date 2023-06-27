%% compute_FSBP_poly
%
% Description: 
%  Function to compute a first- and second-derivative 
%  polynomial SBP operators on Gauss-Lobatto points on [-1,1]
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT:  
%  d :          polynomial degree 
%
% OUTPUT: 
%  D1 : differentiation matrix 
%  D2 : differentiation matrix
%  x :  grid points 
%  P :  diagonal-norm matrix
%  Q :  matrix for boundary correction         

function [D1, D2, x, P, Q] = compute_FSBP_poly( d )

    %% (P1) Compute the first-derivative SBP operator 
    N = d+1; % number of points 
    [x,w,P_aux] = lglnodes(N-1); % Gauss-Lobatto points and weights 
    x = flip(x); % re-order points
    P = diag(w); % diagonal norm matrix 
    D1 = zeros(N,N); % initialize the first-derivative SBP operator 
    for i=1:N
        for j=1:N
            D1(i,j) =dx_lagrange(j,x,x(i)); % assign values to first-derivative operator 
        end
    end
    Q = inv(P)*D1; % get the matrix for boundary correction  

    %% (P2) Get the second-derivative SBP operator 
    D2 = D1*D1; 
    
end