function [solution] = f_max_exp_K3(x)
%{
---------------------------------------------------------------------------
Description:
function handle for the function of Task 1.4 4. 
---------------------------------------------------------------------------
Parameters:
    x: Input vector x (2-dimensional)
---------------------------------------------------------------------------
%}
[N, M] = size(x);
K = 3;

%{
if N ~= 2
    ValueError("Test")
end
%}

for k=1:M
    inner_sum(:, k) = exp(x(1, k)) + exp(x(2, k)) - K;
    solution(k) = max(inner_sum(:, k), 0);
end


end