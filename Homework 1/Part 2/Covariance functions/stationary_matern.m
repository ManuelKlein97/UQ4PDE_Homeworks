function [solution] = stationary_matern(r, nu)
%{
---------------------------------------------------------------------------
Description:
function to return the values of the stationary covariance function
depending on the distance r to the origin 
---------------------------------------------------------------------------
Parameters:
         r:  Discretization scheme of the interval [0,1]
        nu:  Matern Kernel parameter (special cases only)
---------------------------------------------------------------------------
%}

if ~(nu == 0.5 | nu == 1.5 | nu == 2.5 | nu == 'infinity')
    error("Paramter nu has to be either value 0.5, 1.5, 2.5 or string 'infinity'")
end

sigma = sqrt(2);
rho = 0.1;

L = length(r);

solution = zeros([L 1]);

if nu == 0.5
    for l=1:L
        solution(l) = sigma^2*exp(-abs(r(l))/rho);
    end
end

if nu == 1.5
    for l=1:L
        solution(l) = sigma^2*(1 + sqrt(3)*abs(r(l))/rho)*exp(-sqrt(3)*abs(r(l))/rho);
    end
end

if nu == 2.5
    for l=1:L
        solution(l) = sigma^2*(1 + sqrt(5)*abs(r(l))/(rho) + sqrt(3)*abs(r(l))^2/(rho^2))*exp(-sqrt(5)*abs(r(l))/rho);
    end
end

if nu == 'infinity'
    for l=1:L
        solution(l) = sigma^2*exp(-abs(r(l))^2/(2*rho^2));
    end
end



end