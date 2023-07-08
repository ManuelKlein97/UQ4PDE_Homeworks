function [solution] = f_continuous(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the continuous function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (N-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
c = 2.04/N*ones([N 1]);
w = 1/2;

for k=1:M
    inner_sum(:, k) = c.*abs(x(:, k) - w);
    solution(k) = exp(-sum(inner_sum(:, k)));
end
end