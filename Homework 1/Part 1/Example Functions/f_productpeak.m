function [solution] = f_productpeak(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the Product peak function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (N-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
c = 7.25/N*ones([N 1]);
c_inv = 1./c;
w = 1/2*ones([N 1]);

for k=1:M
    inner_prod_vec(:, k) = 1./(c_inv.^2 + (x(:, k) - w).^2);
    solution(k) = prod(inner_prod_vec(:, k));
end


end