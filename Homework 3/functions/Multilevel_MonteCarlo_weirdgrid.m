function [MLMC_estimate, MLMC_variance, Lmax, M_l, workbound] = Multilevel_MonteCarlo(h0, M0, TOL, nu, method, shift)
%{
---------------------------------------------------------------------------
Description:

---------------------------------------------------------------------------
Inputs:
            h0: Stepsize of coarsest level
            M0: Initial number of samples to calculate the variances and
                costs to determine optimal number of samples for each level
           TOL: Tolerance which is supposed to be achieved by the procedure
            nu: random field parameter
        method: QMC or RQMC or (plain) MC
        
Outputs:
      MLMC_estimate: Multilevel Estimation
      MLMC_variance: Multilevel Variance
               Lmax: Maximum number of levels to achieve desired accuracy
                M_l: Optimal number of samples for each level
          workbound:
---------------------------------------------------------------------------
%}

% ------------------------------------------------------------------------
% First we calculate Lmax such that we can assure the Tolerance
% The bias is <=4/3*h^2

%Lmax = ceil(-log(sqrt(3)/2*(TOL/(h0*sqrt(2))))*1/log(2));

Lmax = ceil((-2*log(3*TOL/(8*h0^2))+log(2))/(4*log(2)));

% First we calculate the h-vector containing all the stepsizes for each
% level h(1) = h0, h(l) = h0*2^(-l+1), l>1
% For level L we need discretization of level L+2
h = zeros(1, Lmax+1+2);
for l=1:Lmax+1+2
    h(l) = h0*2^(-l+1);
end


% Now we calculate the finest grid, for that we need the number of points
I = ceil(1./h);
x = linspace(0, 1, I(end)+1);
x_auxiliary = x(2:I(end));


MLMC_estimate = 0;
MLMC_variance = 0;
% ------------------------------------------------------------------------
% Calculate the random field approximation for the finest level Lmax
rho = 0.1;
sigma = sqrt(2);

% Cutoffs for each nu parameter value
if nu == 0.5
    cutoff = 78;
elseif nu == 1.5
    cutoff = 26;
elseif nu == 2.5
    cutoff = 19;
elseif strcmp(nu, 'infinity')
    cutoff = 13;
end
C_finest = zeros(I(end)-1, I(end)-1);
for i=1:I(end)-1
    for j=1:I(end)-1
        if nu == 0.5
            C_finest(i, j) = sigma^2 .* exp( - abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif nu == 1.5
            C_finest(i, j) = sigma^2 .* (1 + sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j))/rho) .* exp( - sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif nu == 2.5
            C_finest(i, j) = sigma^2 .* (1 + sqrt(5)*abs(x_auxiliary(i)-x_auxiliary(j))/rho + sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j)).^2/rho^2) .* exp( - sqrt(5)*abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif strcmp(nu, 'infinity')
            C_finest(i, j) = sigma^2 .* exp(-abs(x_auxiliary(i)-x_auxiliary(j)).^2/(2*rho^2));
        end
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
% Check if cutoff is greater than the number of samples in use
if cutoff > ceil(I(end)/2)
    cutoff = ceil(I(end)/2);
end
% Check if Importance sampling shift is non-zero
if all(shift == 0)
    if nu == 0.5
        cutoff = 29;
        shift = zeros(1, cutoff);
    elseif nu == 1.5
        cutoff = 21;
        shift = zeros(1, cutoff);
    elseif nu == 2.5
        cutoff = 19;
        shift = zeros(1, cutoff);
    elseif strcmp(nu, 'infinity')
        cutoff = 16;
        shift = zeros(1, cutoff);
    end
else
    cutoff = length(shift);
end


% Check method
if strcmp(method, 'MC')
    Y = normrnd(0, 1, [cutoff M0]) + shift';
elseif strcmp(method, 'QMC')
    % Sobolpoint 0 is problematic because norminv(0)=-inf (?)
    sobolpoints = i4_sobol_generate(cutoff, M0+1, 0);
    Y = norminv(sobolpoints(:, 2:M0+1), shift', 1);
elseif strcmp(method, 'RQMC')
    sobolpoints = i4_sobol_generate(cutoff, M0, 0);
    randomized_sobolpoints = sobolpoints + randn(cutoff, M0) + shift';
    uniformed_rand_sp = mod(randomized_sobolpoints, 1);
    Y = norminv(uniformed_rand_sp);
end

% Generate M0 samples of the random field to determine costs and variance
% for optimal number of samples in each level
kappa = zeros(I(end)-1, M0);
randomfield_a = zeros(I(end)-1, M0);
for m=1:M0
    for n=1:cutoff
        kappa(:, m) = kappa(:, m) + sqrt(EVal(n))*Y(n, m)*EVec(:, n);
    end
    randomfield_a(:, m) = exp(kappa(:, m));
end

% ------------------------------------------------------------------------
% Now we will calculate the QoI and get the average sampling 
% costs for each level
% First assemble matrix A, and vector F while collecting the time spend on
% one sample

MC_estimate = zeros(Lmax+1, 1);
variance = zeros(Lmax+1, 1);
costs = zeros(Lmax+1, 1);
f = @(x) 4*pi^2*cos(2*pi*x);
for l=1:Lmax+1
    cost_level = zeros(M0, 1);
    % Array for saving the corresponding indices in x_auxiliary for x_level
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

    QoI = zeros(M0, 1);
    QoI_temp = zeros(M0, 1);
    for m=1:M0
        tic
        A = zeros(I(l), I(l));
        F = zeros(I(l), 1);
        for i=2:I(l)+1
            F(i-1) = f(x_level(i));
        end
        for i=1:I(l)
            if i == 1
                % Case if h/2 has to be used
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
        if l == 1
            % QoI for level 0: sampling g0
            QoI(m) = h(l)*sum(u_h);
        else
            % QoI for level l>0: sampling g_l - g_{l-1}
            QoI_temp(m) = QoI(m);
            QoI(m) = h(l)*sum(u_h) - QoI_temp(m);
        end
        cost_level(m) = toc;
    end
    
    % Calculate MC estimate, variance, avg. cost per sample for QoI
    MC_estimate(l) = 1/M0*sum(QoI);
    variance(l) = 1/(M0-1)*sum((QoI - MC_estimate(l)).^2);
    costs(l) = 1/M0*sum(cost_level);
    %QoI_previous_level = QoI;
end


% ------------------------------------------------------------------------
% Now we can calculate the optimal sample sizes for each level and the
% bound for the optimal work given by slides 436 in the UQ4PDE course
% lecture slides

M_l = zeros(Lmax+1, 1);
for l=1:Lmax+1
    M_l(l) = ceil(TOL^(-2)*sqrt(variance(l)/costs(l))*sum(sqrt(variance.*costs)));
end
workbound = TOL^(-2)*sum(sqrt(variance.*costs))^2 + sum(costs);

% We have to check whether M(l)>M(l+1)
if ~(h0 == 1/2)
    for l=1:Lmax
        if M_l(l) - M_l(l+1)<0
        error("Optimal sample size for each level is not decreasing.")
        end
    end
end

% Print Number of samples if the number of samples of first level exceeds
% threshold (M_0 > 50.000)
if M_l(1) > 5*10^4
    M_l
end
% ------------------------------------------------------------------------
% Now we will calculate a MLMC estimate using the optimal number of samples
% we just calculated

% We can easily implement this procedure using the previous code with
% differnt sample sizes dependend on the level currently worked on

% ------------------------------------------------------------------------
% Treat level 0 (l=1) differently
l = 1;
% Check
if strcmp(method, 'MC')
    Y = normrnd(0, 1, [cutoff M_l(l)]) + shift';
elseif strcmp(method, 'QMC')
    sobolpoints = i4_sobol_generate(cutoff, M_l(l)+1, 0);
    Y = norminv(sobolpoints(:, 2:M_l(l)+1), shift', 1);
elseif strcmp(method, 'RQMC')
    sobolpoints = i4_sobol_generate(cutoff, M_l(l), 0);
    randomized_sobolpoints = sobolpoints + randn(cutoff, M_l(l)) + shift';
    uniformed_rand_sp = mod(randomized_sobolpoints, 1);
    Y = norminv(uniformed_rand_sp);
end
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

% Calculate MC estimate, variance, avg. cost per sample for QoI
MC_est_L0 = 1/M_l(l)*sum(QoI);
MC_var_L0 = 1/(M_l(l)-1)*sum((QoI - MC_est_L0).^2);
time = toc;
disp(['Elapsed time of level 0: ', num2str(time), ' seconds']);

% ------------------------------------------------------------------------
% Now for level 1+ (l=2+)

% We can generate max(M_l) = M_l(2) samples to use in every level
if strcmp(method, 'MC')
    % TODO Maybe add a shift, s.t. we can use Importance sampling
    Y = normrnd(0, 1, [cutoff M_l(l)]) + shift';
elseif strcmp(method, 'QMC')
    sobolpoints = i4_sobol_generate(cutoff, M_l(l)+1, 0);
    Y = norminv(sobolpoints(:, 2:M_l(l)+1), shift', 1);
elseif strcmp(method, 'RQMC')
    sobolpoints = i4_sobol_generate(cutoff, M_l(l), 0);
    randomized_sobolpoints = sobolpoints + randn(cutoff, M_l(l)) + shift';
    uniformed_rand_sp = mod(randomized_sobolpoints, 1);
    Y = norminv(uniformed_rand_sp);
end
kappa = zeros(I(end)-1, M_l(l));
randomfield_a = zeros(I(end)-1, M_l(l));
for m=1:M_l(l)
    for n=1:cutoff
        kappa(:, m) = kappa(:, m) + sqrt(EVal(n))*Y(n, m)*EVec(:, n);
    end
    randomfield_a(:, m) = exp(kappa(:, m));
end

QoI_previous_level = zeros(M_l(l), 1);
for l=1:Lmax+1
    tic
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

    % Calculate MC estimate, variance, avg. cost per sample for QoI
    MC_estimate(l) = 1/M_l(l)*sum(QoI);
    %variance(l) = 1/(M_l(l)-1)*sum((QoI - MC_estimate(l)).^2);
    QoI_previous_level = QoI;
    if l>1
        variance(l-1) = 1/(M_l(l)-1)*sum((QoI - QoI_previous_level - MC_estimate(l)).^2);
    end
    time = toc;
    disp(['Elapsed time of level ', num2str(l), ': ', num2str(time), ' seconds']);
end

% Finally we can write the MLMC estimator
MLMC_estimate = MC_est_L0;
MLMC_variance = MC_var_L0/M_l(1) + sum((1./M_l).*variance);
for l=2:Lmax+1
    MLMC_estimate = MLMC_estimate + (MC_estimate(l) - MC_estimate(l-1));
end
Lmax = Lmax+1;

end