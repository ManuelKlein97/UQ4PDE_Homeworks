function [solution] = a_constant_coeff(x, N, sigma)
%{
---------------------------------------------------------------------------
Description:
function to return a sample of a piecewise constant random field.
---------------------------------------------------------------------------
Parameters:
        x: Function input vector (linspace) or single value
        N: mesh to define random field
    sigma: factor in front of the sum of the field representation
---------------------------------------------------------------------------
%}


% Make sure sigma is chosen such that 1-sigma\sqrt(3)>0 to ensure coercivity
if 1-sigma*sqrt(3) <= 0
    error("Choose sigma such that 1-sigma*sqrt(3)>0, to ensure coercivity.")
end

% Create equidistant grid (N+1 steps) for the indicator function
x_hat = linspace(0, 1, N + 1);

% Sample uniformly on the interval [-sqrt(3), sqrt(3)]

% Set the boundaries of the uniform sampling
a = -sqrt(3);
b = sqrt(3);

% Sample uniformly on [a,b]
Y = a + (b-a)*rand(N, 1);

% Initialize temporary array for indicator function
ind = zeros([length(x) N]);

for l=1:length(x)
    for k=2:N+1
        if (x(l) >= x_hat(k-1) & x(l) <= x_hat(k))
            ind(l, k) = Y(k-1);
        else
            ind(l, k) = 0;
        end
    end
end

for l=1:length(x)
    solution(l) = 1 + sigma*sum(ind(l,:));
end

end