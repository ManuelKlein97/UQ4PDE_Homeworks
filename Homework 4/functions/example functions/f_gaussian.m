function [solution] = f_gaussian(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the Gaussian function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (N-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
c = 7.03/N*ones([N 1]);
w = 1/2;

for k=1:M
    inner_sum(:, k) = (c.^2).*(x(:, k) - w).^2;
    solution(k) = exp(-sum(inner_sum(:, k)));
end

end