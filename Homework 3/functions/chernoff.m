function [P_bound, max_theta] = chernoff(theta, K)
%{
---------------------------------------------------------------------------
Description:
    Calculates the Chernoffbounds from HW3 task 4 for model 2 with
    parameter nu=0.5 and importance sampling shift from HW1.
---------------------------------------------------------------------------
Inputs:
               
Outputs:
     
---------------------------------------------------------------------------
%}
% Setting the parameters for this specific task
% IS shift
shift = [1.77 0 1.56 0 -0.38];
nu = 0.5;
h = [1/4 1/8 1/16 1/32 1/64 1/128 1/256];
I = 1./h;
if K == -5
    M_l = [3748 1440 702 505 344];
elseif K == -10
    M_l = [1362 1274 269 226 115];
elseif K == -20
    M_l = [2782 501 210 170 65];
else
    error("Invalid value for K. Choose -5,-10,-20");
end
cutoff = 5;
I_max = 64;
Ntheta = length(theta);

% ------------------------------------------------------------------------
% First we calculate the eigenvalues and eigenvectors for our randomfield
% representation
x = linspace(0, 1, I(end)+1);
x_auxiliary = x(2:I(end));
rho = 0.1;
sigma = sqrt(2);
C_finest = zeros(I(end)-1, I(end)-1);
for i=1:I(end)-1
    for j=1:I(end)-1
        C_finest(i, j) = sigma^2 .* exp( - abs(x_auxiliary(i)-x_auxiliary(j))/rho);
    end
end
C_finest = h(end).*C_finest;
[unsorted_EVec, unsorted_EVal] = eig(C_finest);
for i = 1:I(end)-1
    % Norm the EV
    unsorted_EVec(:, i) = unsorted_EVec(:, i)./(sqrt(h(end))*norm(unsorted_EVec(:, i)));
end
% Sort the Eigenvalues/Eigenvectors
EVec = zeros(I(end)-1, I(end)-1); 
EVal = zeros(I(end)-1, 1);
for m = 1:I(end)-1
    EVec(:, m) = unsorted_EVec(:, I(end)-m);
    EVal(m) = unsorted_EVal(I(end)-m, I(end)-m);
end


f = @(x) 4*pi^2*cos(2*pi*x);
% Now we can sample QoI for level 0
% First the random field
l = 1;
Y = normrnd(0, 1, [cutoff M_l(l)]) + shift';
kappa = zeros(I(end)-1, M_l(l));
randomfield_a = zeros(I(end)-1, M_l(l));
for m=1:M_l(l)
    for n=1:cutoff
        kappa(:, m) = kappa(:, m) + sqrt(EVal(n))*Y(n, m)*EVec(:, n);
    end
    randomfield_a(:, m) = exp(kappa(:, m));
end

indices = zeros(I(l), 1);
x_level = zeros(I(l), 1);
start = find(x_auxiliary == h(l)/2);
x_level(1) = x_auxiliary(start);
indices(1) = start;
for j=2:I(l)
    index = find(x_auxiliary == x_level(j-1) + h(l));
    indices(j) = index;
    x_level(j) = x_auxiliary(index);
end
x_level = [0; x_level; 1];
x_A = zeros(length(x_level)-1, 1);
A_index = zeros(length(x_level)-1, 1);
for i=1:length(x_level)-1
    x_A(i) = (x_level(i) + x_level(i+1))/2;
    A_index(i) = find(x_auxiliary == x_A(i));
end
QoI = zeros(M_l(l), 1);
for m=1:M_l(l)
    A = zeros(I(l), I(l));
    F = zeros(I(l), 1);
    for i=2:I(l)+1
        F(i-1) = f(x_level(i));
    end
    for i=1:I(l)
        if i == 1
            A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/(h(l)/2)^2;
            A(i, i+1) = -randomfield_a(A_index(i+1), m)/(h(l)/2)^2;
        elseif i == I(l)
            A(i, i-1) = - randomfield_a(A_index(i), m)/(h(l)/2)^2;
            A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/(h(l)/2)^2;
        else
            A(i, i-1) = - randomfield_a(A_index(i), m)/h(l)^2;
            A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/h(l)^2;
            A(i, i+1) = -randomfield_a(A_index(i+1), m)/h(l)^2;
        end
    end
    u_h = A\F;
    % QoI for level 0: sampling g0
    QoI(m) = h(l)*sum(u_h);
end

% First part of exp_estimate for each theta value, needs to be updated (add
% other levels differences)
exp_est = zeros(Ntheta, 1);
for n=1:Ntheta
    % This needs to be MLMC estimator instead of plain MC
    exp_est(n) = 1/M_l(l)*sum(exp(-theta(n)*QoI));
end



% Generate M_l(1) random field samples with stepsize h_Lmax
Y = normrnd(0, 1, [cutoff M_l(l)]) + shift';
kappa = zeros(I(end)-1, M_l(l));
randomfield_a = zeros(I(end)-1, M_l(l));
for m=1:M_l(l)
    for n=1:cutoff
        kappa(:, m) = kappa(:, m) + sqrt(EVal(n))*Y(n, m)*EVec(:, n);
    end
    randomfield_a(:, m) = exp(kappa(:, m));
end

% Now calculate QoIs for each level
Lmax = 5; % Lmax is theoretically 4 but we start indexing at 1 instead of 0
for l=1:Lmax
    indices = zeros(I(l), 1);
    x_level = zeros(I(l), 1);
    start = find(x_auxiliary == h(l)/2);
    x_level(1) = x_auxiliary(start);
    indices(1) = start;
    for j=2:I(l)
        index = find(x_auxiliary == x_level(j-1) + h(l));
        indices(j) = index;
        x_level(j) = x_auxiliary(index);
    end
    x_level = [0; x_level; 1];
    x_A = zeros(length(x_level)-1, 1);
    A_index = zeros(length(x_level)-1, 1);
    for i=1:length(x_level)-1
        x_A(i) = (x_level(i) + x_level(i+1))/2;
        A_index(i) = find(x_auxiliary == x_A(i));
    end
    QoI = zeros(M_l(l), 1);
    for m=1:M_l(l)
        A = zeros(I(l), I(l));
        F = zeros(I(l), 1);
        for i=2:I(l)+1
            F(i-1) = f(x_level(i));
        end
        for i=1:I(l)
            if i == 1
                A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/(h(l)/2)^2;
                A(i, i+1) = -randomfield_a(A_index(i+1), m)/(h(l)/2)^2;
            elseif i == I(l)
                A(i, i-1) = - randomfield_a(A_index(i), m)/(h(l)/2)^2;
                A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/(h(l)/2)^2;
            else
                A(i, i-1) = - randomfield_a(A_index(i), m)/h(l)^2;
                A(i, i) = (randomfield_a(A_index(i), m) + randomfield_a(A_index(i+1), m))/h(l)^2;
                A(i, i+1) = -randomfield_a(A_index(i+1), m)/h(l)^2;
            end
        end
        u_h = A\F;
        QoI(m) = h(l)*sum(u_h);
    end
    if l>1   % If l>1 we can do the difference operator evaluation
        for n=1:Ntheta
            % This needs to be MLMC estimator instead of plain MC
            exp_est(n) = exp_est(n) + 1/M_l(l)*sum(exp(-theta(n)*QoI)) - 1/M_l(l-1)*sum(exp(-theta(n)*QoI_prev));
        end
    end
    QoI_prev = QoI;
end
I_est = zeros(Ntheta, 1);
for n=1:Ntheta   
    Lambda(n) = log(exp_est(n));
    I_est(n) = theta(n)*(-K) - log(exp_est(n));
end

[I_est_sup, max_theta] = max(I_est);

max_theta = theta(max_theta);

P_bound = exp(-I_est_sup);

end