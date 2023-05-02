function [solution] = targetfunction_K3(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the target function of the shift dilation optimization 
problem given in Homework 1.
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (2-dimensional)
---------------------------------------------------------------------------
%}

[N, M] = size(x);

for k=1:M
    inner(:, k) = max(exp(x(1, k)) + exp(x(2, k)) - 3, 0);
    solution(k) = inner(:, k) *  (1/(2*pi)) * exp(-0.5*x(1, k).^2) * exp(-0.5*x(2, k).^2);
end


end