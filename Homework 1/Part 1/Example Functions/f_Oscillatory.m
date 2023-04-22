function [solution] = f_Oscillatory(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the Oscillator function f given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (N-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
c = 9/N*ones([1 N]);
w_1 = 1/2;

for i=1:M
    solution(i) = cos(2*pi*w_1 + dot(c, x(:, i)));
end

end