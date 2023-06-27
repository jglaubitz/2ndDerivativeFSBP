%% dx_lagrange  
%
% Description: 
%  Function to compute the derivative of the Lagrange basis polynomials  
%
% Author: Jan Glaubitz 
% Date: June 27, 2023
% 
% INPUT: 
%  j :          index of the Lagrange basis polynomial  
%  points :     interpolation points 
%  z :          argument for which derivative is evaluated 
%
% OUTPUT: 
%  y :	        nodal value of the derivative at z

function y = dx_lagrange(j,points,z)

    y = 0;
    n = length(points);
    for l=1:n
        if not(l==j)
            k = 1/(points(j)-points(l));
            for m=1:n
                if not(m==j) && not(m==l)
                    k = k*(z-points(m))/(points(j)-points(m));
                end
            end
            y = y + k;
        end
    end

end