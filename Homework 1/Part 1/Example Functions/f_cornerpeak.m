function [solution] = f_cornerpeak(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the corner peak function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (N-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
c = 1.85/N*ones([N 1]);


for k=1:M
    inner_prod_vec(:, k) = 1 + dot(c, x(:, k));
    solution(k) = inner_prod_vec(:, k)^(-N+1);
end


end