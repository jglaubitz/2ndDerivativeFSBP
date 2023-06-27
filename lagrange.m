%% lagrange  
%
% Description: 
%  Function to compute and evaluate the Lagrange basis polynomials  
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT: 
%  j :          index of the Lagrange basis polynomial  
%  points :     interpolation points 
%  z :          argument for which Lagrange basis polynomials are evaluated 
%
% OUTPUT: 
%  y :	        nodal value of the function values at z

function y = lagrange(j,points,z)

    y = 1;
    n = length(points);
    for l=1:n
        if not(l==j) 
            y = y*(z-points(l))/(points(j)-points(l));
        end
    end

end