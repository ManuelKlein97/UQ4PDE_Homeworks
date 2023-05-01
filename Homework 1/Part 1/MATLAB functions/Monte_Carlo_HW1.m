function [MC_estimate, MC_var] = Monte_Carlo_HW1(f, N, M);
%{
---------------------------------------------------------------------------
Description:
This function takes an input which is a function handle and computes the
Monte Carlo estimate of the function f in the hypercube [0,1]^N. f is
supposed to be a n-dimensional function, i.e. its input is an array
x=(x_1,...,x_N).
---------------------------------------------------------------------------
Parameters:
    f: function handle
    N: dimension of the function handle inputs
    M: number of samples used in MC estimate
---------------------------------------------------------------------------
%}
samples = rand(N, M);
y_values = f(samples);
MC_estimate = 1/M*sum(y_values);

MC_var = (1/(M-1))*sum((y_values - MC_estimate).^2);

end