function [MC_est, MC_var] = Monte_Carlo_HW1_Gaussian_Samples_shifted(f, M, mu1, mu2)
%{
---------------------------------------------------------------------------
Description:
This function takes an input which is a function handle and computes the
Monte Carlo estimate of the function f using std. Gaussian samples. f is
supposed to be a 2-dimensional function, i.e. its input is an array
x=(x_1,x_2).
---------------------------------------------------------------------------
Parameters:
    f: function handle
    M: number of samples used in MC estimate
    mu1: Shift coming from the shift dilation technique in IS
    mu2: Shift coming from the shift dilation technique in IS
---------------------------------------------------------------------------
%}
sample1 = normrnd(mu1, 1, [1, M]);
sample2 = normrnd(mu2, 1, [1, M]);
samples = [sample1; sample2];
y_values = f(samples);
MC_est = (1/M)*sum(y_values);

MC_var = (1/(M-1))*sum((y_values - MC_est).^2);

end