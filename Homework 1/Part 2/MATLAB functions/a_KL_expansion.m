function [solution] = a_KL_expansion(x, nu, N)
%{
---------------------------------------------------------------------------
Description:
function to generate a sample path of the random field a(x,w) of model 2
(Log normal) using a Karhunen-Loeve expansion on a given grid x, nu
---------------------------------------------------------------------------
Parameters:
         x: Discretization scheme of the interval [0,1] (size I+1)
        nu: Matern Kernel parameter (special cases only)
         N: N-term truncation
---------------------------------------------------------------------------
%}

% Generation of grid {x_{i+1/2}}_i

I = length(x) - 1;
x_half = zeros([I+1 1]);
for k=1:I
    x_half(k) = (x(k+1) + x(k))/2; 
end

% Generate the covariance matrix
C = Matern_specialcases_only(x_half, nu);

% Solve Eigenvalue problem
[EVec, EVal] = eig(C);

SV = zeros([I 1]); % Singular values array
LinfNorm_EVec = zeros([I 1]); % Inf-Norm array of EVec's

% Calculation of singular values & rearangement of order
EVec = flip(EVec, 2);
for k=1:I
    SV(k) = sqrt(EVal(I+1-k, I+1-k));
    LinfNorm_EVec(k) = norm(EVec(:, k), "inf");
end

% Generate random samples from std Gaussian
Y = normrnd(0, 1, [N 1]);

% Generate kappa array
inner = zeros([I+1 N]);
for k=1:N
    inner(:, k) = SV(k)*Y(k).*EVec(:, k);
end

kappa = sum(inner, 2);
solution = exp(kappa);

end