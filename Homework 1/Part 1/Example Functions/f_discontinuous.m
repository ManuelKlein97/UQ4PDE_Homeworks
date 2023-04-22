function [solution] = f_discontinuous(x)
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
c = 4.3/N*ones([N 1]);
w = [pi/4 pi/5];

for k=1:M
    if ((x(1, k) > w(1)) | (x(2, k)>w(2)))
        solution(k) = 0;
    else
        solution(k) = exp(dot(c, x(:, k)));
    end
end

end