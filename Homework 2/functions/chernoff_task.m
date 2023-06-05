function [P_bound, max_theta] = chernoff_task(nu, M, theta, I, X)

h = 1./I;

Ntheta = length(theta);
NX = length(X);

if nu == 0.5
    cutoff = 29;
elseif nu == 1.5
    cutoff = 21;
elseif nu == 2.5
    cutoff = 19;
elseif strcmp(nu, 'infinity')
    cutoff = 16;
end

% Generate normal samples (2 times)
sobolpoints1 = i4_sobol_generate(cutoff, M, 0);
rng(1); % Seed for reproducability
randomized_sobolpoints1 = sobolpoints1 + randn(cutoff, M);
uniformed_rand_sp1 = mod(randomized_sobolpoints1, 1);
y = norminv(uniformed_rand_sp1);


sobolpoints2 = i4_sobol_generate(cutoff, M, 0);
rng(2); % Seed for reproducability
randomized_sobolpoints2 = sobolpoints2 + randn(cutoff, M);
uniformed_rand_sp2 = mod(randomized_sobolpoints2, 1);
z = norminv(uniformed_rand_sp2);

% Now we can represent the random field a using the Fourier series
if strcmp(nu, 'infinity')
    nu = 0;
end
fourier = fourierkoeff(nu, cutoff);

% Construction of F for normal mesh size
F = zeros(I-1, 1);
for i=1:I-1
    F(i) = f(i*h);
end

QoI = zeros(1, M);
exp_estimate = zeros(Ntheta, 1);
Lambda = zeros(Ntheta, 1);
I_est = zeros(NX, Ntheta);
for m=1:M
    % Construction of A for normal mesh size
    A = assemble_A(@a2_3, I, h, fourier, cutoff, y(:, m), z(:, m));
    uh = A\F;
    QoI(m) = h*sum(uh);
end
for n=1:Ntheta   
    exp_estimate(n) = 1/M*sum(exp(-theta(n)*QoI(1:M)));
    Lambda(n) = log(exp_estimate(n));
    for i=1:NX
        I_est(i, n) = theta(n)*X(i) - Lambda(n);
    end
end

I_est_sup = zeros(NX, 1);
max_theta = zeros(NX, 1);
% Now optimizing for the supremum
for i=1:NX
    [I_est_sup(i), max_theta(i)] = max(I_est(i, :));
end
max_theta = theta(max_theta);

P_bound = exp(-I_est_sup);

end