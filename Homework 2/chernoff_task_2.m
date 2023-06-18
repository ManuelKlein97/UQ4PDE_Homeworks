function [P_bound, max_theta] = chernoff_task_2(nu, I, K, theta, M)

Ntheta = length(theta);
theta = 0.4712*ones(Ntheta, 1);

x = linspace(0, 1, I+1);
x_half = zeros([I+1 1]);
for k=1:I
    x_half(k) = (x(k+1) + x(k))/2; 
end

Cov = Matern_specialcases_only(x_half, nu);
[EVec, EVal] = eig(Cov);
SV = zeros([I 1]);
b = zeros(size(EVec));
max_norm_b = zeros([length(EVec) 1]);
EVec = flip(EVec, 2);
for k=1:I
    SV(k) = sqrt(EVal(I+1-k, I+1-k));
    b(:, k) = EVec(:, k) * SV(k);
    max_norm_b(k) = max(b(:, k));
end
b0 = min(b(:, 1));
q0 = (0.5*sqrt(8)*pi^2)/b0;
exp_value_abs_std_gauss = sqrt(2/pi);
innersum = sum(max_norm_b)*exp_value_abs_std_gauss;
innerI = log(K/q0)-innersum;


exp_estimate = zeros(Ntheta, 1);
Lambda = zeros(Ntheta, 1);
I_est = zeros(Ntheta, 1);
Y = randn(I);
for k=1:Ntheta
    for n=1:I
        exp_estimate(k) = exp_estimate(k) + log(1/M*sum(exp(-theta(k)*max_norm_b(n)*(abs(Y(n)) - sqrt(2/pi)))));
    end
    Lambda(k) = exp_estimate(k);
    
    I_est(k) = theta(k)*innerI - Lambda(k);
    
end

% Now optimizing for the supremum
[I_est_sup, max_theta] = max(I_est);

max_theta = theta(max_theta);

P_bound = exp(-I_est_sup);



end