function [P_bound, opt_theta] = chernoff_prod(theta, K, nu)
%{
---------------------------------------------------------------------------
Description:
    Calculates the Chernoffbounds from HW3 task 4 for model 2 with
    parameter nu=0.5 for my suggestion of the approximation of the
    exponential moments using products instead of sums.
---------------------------------------------------------------------------
Inputs:
               
Outputs:
     
---------------------------------------------------------------------------
%}
h = [1/4 1/8 1/16 1/32 1/64 1/128 1/256];
I = 1./h;

if nu == 0.5
    %cutoff = 29;
    if K == -5
        M_l = 20000*ones(5, 1);
        %M_l = [3748 1440 702 505 344];
    elseif K == -10

        M_l = 20000*ones(5, 1);
        %M_l = [1362 1274 269 226 115];
    elseif K == -20
        M_l = 20000*ones(5, 1);
        %M_l = [2782 501 210 170 65];
    else
        error("Invalid value for K. Choose -5,-10,-20");
    end
elseif nu == 0
    %cutoff = 16;
    % Maybe change M_l values here
    if K == -5
        M_l = [3748 1440 702 505 344];
    elseif K == -10
        M_l = [1362 1274 269 226 115];
    elseif K == -20
        M_l = [2782 501 210 170 65];
    else
        error("Invalid value for K. Choose -5,-10,-20");
    end
else
    error("Invalid nu value.")
end

LMax = length(M_l);
Ntheta = length(theta);

f = @(x) 4*pi^2*cos(2*pi*x);

g = zeros(LMax, Ntheta); % E[exp(-theta(Ntheta)*QoI(LMax))]

% -----------------------------------------------------------------------
% set cutoff to number of discretization points whenever I(l)+1<cutoff
cutoff = I(1);
% Level 0
Y = normrnd(0, 1, [cutoff M_l(1)]);
[EVal, EVec] = getEigen(nu, cutoff, I(1));
randomfield_a = Karhunen_Loeve_N(Y, EVal, EVec, cutoff);
F = zeros(cutoff, 1);
QoI = zeros(M_l(1), 1);
[x_level, ~] = grid_level(h(1));
for i=1:length(x_level)
    F(i) = f(x_level(i));
end
for m=1:M_l(1)
    A = getA(randomfield_a(:, m));
    u_h = A\F;
    QoI(m) = h(1)*sum(u_h);
end
% Generate E[exp(-theta*Q(u_{h_0}))]
for n=1:Ntheta
    g(1, n) = 1/M_l(1) * sum(exp(-theta(n)*QoI));
end

% Now level 1+

for l=2:LMax
    cutoff1 = I(l-1);
    cutoff2 = I(l);
    if nu == 0.5
        if cutoff1 > 29
            cutoff1 = 29;
        end
        if cutoff2 > 29
            cutoff2 = 29;
        end
    end
    if nu == 0
        if cutoff1 > 29
            cutoff1 = 29;
        end
        if cutoff2 > 29
            cutoff2 = 29;
        end
    end
    Y = normrnd(0, 1, [cutoff2 M_l(l)]);
    [EVal1, EVec1] = getEigen(nu, cutoff1, I(l-1));
    [EVal2, EVec2] = getEigen(nu, cutoff2, I(l));
    randomfield_a1 = Karhunen_Loeve_N(Y, EVal1, EVec1, cutoff1);
    randomfield_a2 = Karhunen_Loeve_N(Y, EVal2, EVec2, cutoff2);
    F1 = zeros(cutoff1, 1);
    F2 = zeros(cutoff2, 1);
    QoI1 = zeros(M_l(l), 1);
    QoI2 = zeros(M_l(l), 1);
    [x_level1, ~] = grid_level(h(l-1));
    [x_level2, ~] = grid_level(h(l));
    for i=1:length(x_level1)
        F1(i) = f(x_level1(i));
    end
    for i=1:length(x_level2)
        F2(i) = f(x_level2(i));
    end
    for m=1:M_l(l)
        A1 = getA(randomfield_a1(:, m));
        u_h1 = A1\F1;
        A2 = getA(randomfield_a2(:, m));
        u_h2 = A2\F2;
        QoI1(m) = h(l-1)*sum(u_h1);
        QoI2(m) = h(l)*sum(u_h2);
    end
    % Calculate expected value of 
    for n=1:Ntheta
        g(l, n) = - 1/M_l(l) * sum(exp(-theta(n)*QoI1)) + 1/M_l(l) * sum(exp(-theta(n)*QoI2));
    end
end

% Using product representation of E[exp(-theta*QoI)] we get
Lambda = log(sum(g, 1));
I_sup = theta*(-K) - Lambda;
[I_sup_optimaltheta, theta_index] = max(I_sup);
opt_theta = theta(theta_index);

% Finally our bound:
P_bound = exp(-I_sup_optimaltheta);

end