function [coeff] = getLegendreCoefficients1D(f, w)
%{
---------------------------------------------------------------------------
Description:
Function to return the legendre coefficients for the truncated legendre
series of the function f with truncation parameter w
---------------------------------------------------------------------------
Parameters:
    f: f function-handle
    w: truncation parameter
---------------------------------------------------------------------------
%}

% Account for p=0 and matlabs useless 1-start
w = w + 1;

% Calculate the p-th coeffiecient (p <= w)
coeff = zeros(1, w);

for p=1:w
    integrand = @(x) sqrt(2*(p-1)+1)*legendreP(p-1, 2*x-1).*f(x);
    coeff(p) = integral(integrand, 0, 1);
end

end