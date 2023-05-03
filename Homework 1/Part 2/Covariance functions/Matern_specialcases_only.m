function [solution] = Matern_specialcases_only(x, nu)
%{
---------------------------------------------------------------------------
Description:
function to return the covariance matrix of the discretized
---------------------------------------------------------------------------
Parameters:
         x:  Discretization scheme of the interval [0,1]
        nu:  Matern Kernel parameter (special cases only)
---------------------------------------------------------------------------
%}

if ~(nu == 0.5 | nu == 1.5 | nu == 2.5 | nu == 'infinity')
    error("Paramter nu has to be either value 0.5, 1.5, 2.5 or string 'infinity'")
end

% Set values for the other parameters (sigma^2=2, rho=0.1)
sigma = sqrt(2);
rho = 0.1;

solution = zeros([length(x) length(x)]);

if nu == 0.5
    for l=1:length(x)
        for k=1:length(x)
            solution(l, k) = sigma^2*exp(-abs(x(l) - x(k))/rho);
        end
    end
end

if nu == 1.5
    for l=1:length(x)
        for k=1:length(x)
            solution(l, k) = sigma^2*(1 + sqrt(3)*abs(x(l) - x(k))/rho)*exp(-sqrt(3)*abs(x(l) - x(k))/rho);
        end
    end
end

if nu == 2.5
    for l=1:length(x)
        for k=1:length(x)
            solution(l, k) = sigma^2*(1 + sqrt(5)*abs(x(l) - x(k))/(rho))*(1 + sqrt(3)*abs(x(l) - x(k))/(rho^2))*exp(-sqrt(3)*abs(x(l) - x(k))/rho);
        end
    end
end

if nu == 'infinity'
    for l=1:length(x)
        for k=1:length(x)
            solution(l, k) = sigma^2*exp(-abs(x(l) - x(k))^2/(2*rho^2));
        end
    end
end


end