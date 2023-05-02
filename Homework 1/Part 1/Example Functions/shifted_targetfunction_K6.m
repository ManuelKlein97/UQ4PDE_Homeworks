function [solution] = shifted_targetfunction_K6(x, mu1, mu2)
%{
---------------------------------------------------------------------------
Description:
function handle for the shifted function coming from the importance
sampling shift dilation technique
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (2-dimensional) (Samples)
    mu1, mu2: Shifts
---------------------------------------------------------------------------
%}
[N, M] = size(x);
front = f_max_exp_K6(x);
for k=1:M
    back_1(k) = exp(-1/2*(x(1, k).^2 + x(2, k).^2));
    back_2(k) = exp(-1/2*((x(1, k)- mu1).^2 + (x(2, k)- mu2).^2));
    solution(k) = front(k)*(back_1(k)/back_2(k));
end
end