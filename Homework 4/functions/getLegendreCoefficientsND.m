function [coeff, indxset, Dprod] = getLegendreCoefficientsND(f, w, N, type, M)
%{
---------------------------------------------------------------------------
Description:
Function to return the legendre coefficients for the truncated legendre
series of the function f with truncation parameter w in N dimensions
---------------------------------------------------------------------------
Parameters:
    f: f function-handle
    w: truncation parameter
    N: inputdimension of function f
 type: Polynomial space type ('TP', 'HC', 'TD')
    M: Number of random samples used

Output:
    coeff: Legendre coefficients for multidimensional representation
    indxset: corresponding multi-index set
    Dprod: D^T.D to check stability and condition
---------------------------------------------------------------------------
%}
indxset = generateMultiIndexSet(N, w, type);
[indxsize, ~] = size(indxset);
rng(1)
y = rand([N M]);
phi = zeros([M 1]);
for i=1:M
    phi(i) = f(y(:, i));
end

D = ones(M, indxsize);
for p=1:indxsize
    for i=1:M
        for n=1:N
            D(i, p) = D(i, p)*legendreP(indxset(p, n), 2*y(n, i)-1);
        end
    end
end
Dtrans = transpose(D);
Dprod = mtimes(Dtrans, D);
Dinv = inv(Dprod);
Dfinal = mtimes(Dinv, Dtrans);
coeff = mtimes(Dfinal, phi);
end