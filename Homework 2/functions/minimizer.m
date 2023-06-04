function mu = minimizer(start, K, j, I, h, nu)
    [EWsort, EVsort] = eigC(I, h, nu);
    rho = 0.1;
    c_alpha = 1.96;
    sigma = sqrt(2);
    lb = [];
    ub = [];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    c = [];
    fminimizer = @(x) standardNormal(x)*(-1);
    options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
    mu = fmincon(fminimizer, start, A, b, Aeq, beq, lb, ub, @(x) nzcon(rho, sigma, c_alpha, I, h, nu, x, K, EWsort, EVsort))
    %v=-CalculateQoIforIS(c_alpha, I, h, nu, 5, mu, K, EWsort, EVsort);
end
