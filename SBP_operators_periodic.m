
function [ D1, D2, H ] = SBP_operators_periodic( n, dx, order )

% SBP operators of orders 2, 4, 6 and 8 for the first and second derivatives.

% Inputs:
% n - number of spatial grid pts
% dx - step size
% order - order of accuracy (only for 2,4,6,8)

% Outputs:
% D1 - First derivative operator
% D2 - Second derivative operator (D = P^{-1}M)
% H - The norm operator (denoted P in some papers)


e = ones(n,1);

if order==2
    
    vec = [ 0, 1/2, zeros(1,n-3), -1/2 ];
    D1 = 1/dx*circulant(vec,1); 
    
    %%%%%%
    vec = [ -2, 1, zeros(1,n-3), 1 ];
    D2 = 1/(dx^2)*circulant(vec,1);
    
    H = dx*spdiags([e],0,n,n);
    H(1,1) = dx*1/2;
    H(n,n) = dx*1/2;
  
    
elseif order==4
    
    vec = [ 0, 2/3, -1/12, zeros(1,n-5), 1/12, -2/3 ];
    D1 = 1/dx*circulant(vec,1);
    
    %%%%%%
    vec = [ -5/2, 4/3, -1/12, zeros(1,n-5), -1/12, 4/3 ];
    D2 = 1/(dx^2)*circulant(vec,1);
    
    H = dx*spdiags(e,0,n,n);
    H(1,1) = dx*17/48;
    H(2,2) = dx*59/48;
    H(3,3) = dx*43/48;
    H(4,4) = dx*49/48;
    H(n,n) = H(1,1);
    H(n-1,n-1) = H(2,2);
    H(n-2,n-2) = H(3,3);
    H(n-3,n-3) = H(4,4);
    
    
elseif order==6
    
    vec = [ 0, 3/4, -3/20, 1/60, zeros(1,n-7), -1/60, 3/20, -3/4 ];
    D1 = 1/dx*circulant(vec,1);
    
    %%%%%%
    vec = [ -49/18, 3/2, -3/20, 1/90, zeros(1,n-7), 1/90, -3/20, 3/2 ];
    D2 = 1/(dx^2)*circulant(vec,1);
    
    %%%%%%
    H = dx*spdiags([e],0,n,n); 
    H(1,1) = dx*13649/43200;
    H(2,2) = dx*12013/8640;
    H(3,3) = dx*2711/4320;
    H(4,4) = dx*5359/4320;
    H(5,5) = dx*7877/8640;
    H(6,6) = dx*43801/43200;
    H(n,n) = H(1,1);
    H(n-1,n-1) = H(2,2);
    H(n-2,n-2) = H(3,3);
    H(n-3,n-3) = H(4,4);
    H(n-4,n-4) = H(5,5);
    H(n-5,n-5) = H(6,6);
 
    
else
    disp('Only order 2, 4, and 6 implemented here.')
end