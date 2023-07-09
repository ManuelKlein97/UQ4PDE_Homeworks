function [coeff] = getLegendreCoefficientsND(f, w, N, type, M)
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
---------------------------------------------------------------------------
%}
indxset = generateMultiIndexSet(N, w, type);
y = rand([N M])
scatter(y(1,:), y(2,:))
y = sort(y)
figure()
scatter(y(1,:), y(2,:))
phi = zeros([1 M]);
for i=1:M
    phi(i) = f(y(:, i));
end
phi



end