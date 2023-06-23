function [MLMC_estimate, MLMC_variance] = Multilevel_MC(h0, LMax, M_l, nu)
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
Outputs:
     MLMC_estimate:
     MLMC_variance:             
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
cutoff = I(1)+1;
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
[g(1), V(1)] = PlainMC(QoI);

% Now level 1+

for l=2:LMax+1
    % For level 3 + and nu=0.5 set cutoff 
    if strcmp(nu, 'infinity') && l>3
        cutoff = 13;
    end
    if nu == 0.5 && l>4
        cutoff = 29;
    end
    % Now the level l MC (difference)
    Y = normrnd(0, 1, [cutoff M_l(l)]);
    [EVal1, EVec1] = getEigen(nu, cutoff, I(l-1));
    [EVal2, EVec2] = getEigen(nu, cutoff, I(l));
    randomfield_a1 = Karhunen_Loeve_N(Y, EVal1, EVec1, cutoff);
    randomfield_a2 = Karhunen_Loeve_N(Y, EVal2, EVec2, cutoff);
    [~, x_level1] = grid_level(h(l-1));
    %x_level1 = linspace(0, 1, 1/h(l-1) + 1);
    [~, x_level2] = grid_level(h(l));
    %x_level2 = linspace(0, 1, 1/h(l) + 1);
    F1 = zeros(length(x_level1), 1);
    F2 = zeros(length(x_level2), 1);
    QoI = zeros(M_l(l), 1);
    for i=1:length(x_level1)
        F1(i) = f(x_level1(i));
    end
    for i=1:length(x_level2)
        F2(i) = f(x_level2(i));
    end
    for m=1:M_l(l)
        A1 = getA(randomfield_a1(:, m));
        A2 = getA(randomfield_a2(:, m));
        u_h1 = A1\F1;
        u_h2 = A2\F2;
        QoI(m) = h(l-1)*sum(u_h1) - h(l)*sum(u_h2);
    end
    [g(l), V(l)] = PlainMC(QoI);
end

% Finally we can calculate the MLMC estimate and MLMC variance
MLMC_estimate = sum(g);
MLMC_variance = sum((1./M_l').*V);
end