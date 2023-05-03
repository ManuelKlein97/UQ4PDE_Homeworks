function [solution] = a_matern(x, nu)
%{
---------------------------------------------------------------------------
Description:
function handle for the Oscillator function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
         x: Discretization scheme of the interval [0,1]
        nu: Matern Kernel parameter (special cases only)
---------------------------------------------------------------------------
%}

% Sample the covariance function using the discr. scheme
C = Matern_specialcases_only(x, nu);

% Compute the Cholesky decomposition of C (for nu=0.5, 1.5)
if (nu == 0.5 | nu == 1.5)
    A = chol(C,'lower');
end

if (nu == 2.5 | nu == 'infinity')
    [U,S,V] = svd(C);
    S_sqrt = sqrt(S);
    A = U*S_sqrt;
end

% Sample a normal RV
Z = normrnd(0, 1, [length(x) 1]);

% Get final sample by transformation (since mean zero)
solution = exp(mtimes(A, Z));

end