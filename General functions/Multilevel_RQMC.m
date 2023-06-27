function [MLMC_estimate, MLMC_variance, level_var] = Multilevel_RQMC(h0, LMax, M_l, nu, s)
%{
---------------------------------------------------------------------------
Description:
    Returns Multilevel Monte Carlo estimate and variance.
---------------------------------------------------------------------------
Inputs:
          h0: Initial stepsize
        LMax: Maximum number of levels
         M_l: Samplesize vector [optimally (Lmax x 1)]
          nu: Model 2 parameter
           s: Number of random shifts
Outputs:

    sobolpoints = i4_sobol_generate(cutoff, M_l(1), 0);
    randomized_sobolpoints = sobolpoints + randn(cutoff, M_l(1));
    uniformed_rand_sp = mod(randomized_sobolpoints, 1);
    Y = norminv(uniformed_rand_sp);

     MLMC_estimate:
     MLMC_variance:       
             costs: costs per level per sample
---------------------------------------------------------------------------
%}
h = zeros(LMax+1, 1);
h(1) = h0;
for l=2:LMax+1
    h(l) = h0*2^(-l+1);
end
I = 1./h;
f = @(x) 4*pi^2*cos(2*pi*x);

g = zeros(LMax+1, 1);
V = zeros(LMax+1, 1);

% ------------------------------------------------------------------------
cutoff = I(1);
% Level 0

sobolpoints = i4_sobol_generate(cutoff, M_l(1), 0);
randomized_sobolpoints = sobolpoints + randn(cutoff, M_l(1));
uniformed_rand_sp = mod(randomized_sobolpoints, 1);
Y = norminv(uniformed_rand_sp);


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
[g(1), V(1)] = PlainMC(QoI);

% Now level 1+

for l=2:LMax+1
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

    sobolpoints = i4_sobol_generate(cutoff2, M_l(l), 0);
    randomized_sobolpoints = sobolpoints + randn(cutoff2, M_l(l));
    uniformed_rand_sp = mod(randomized_sobolpoints, 1);
    Y = norminv(uniformed_rand_sp);

    [EVal1, EVec1] = getEigen(nu, cutoff1, I(l-1));
    [EVal2, EVec2] = getEigen(nu, cutoff2, I(l));
    randomfield_a1 = Karhunen_Loeve_N(Y, EVal1, EVec1, cutoff1);
    randomfield_a2 = Karhunen_Loeve_N(Y, EVal2, EVec2, cutoff2);
    F1 = zeros(cutoff1, 1);
    F2 = zeros(cutoff2, 1);
    QoI = zeros(M_l(l), 1);
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
        QoI(m) = h(l)*sum(u_h2) - h(l-1)*sum(u_h1);
    end
    [g(l), V(l)] = PlainMC(QoI);
end
g
level_var = V;
MLMC_estimate = sum(g);
MLMC_variance = sum(1./M_l.*V);

end