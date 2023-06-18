function [MLMC_estimate, MLMC_variance, Lmax, M_l] = Multilevel_MonteCarlo(h0, M0, TOL, nu, method)
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
---------------------------------------------------------------------------
%}

% ------------------------------------------------------------------------
% First we calculate Lmax such that we can assure the Tolerance
% The bias is <=4/3*h^2
Lmax = ceil(-log(sqrt(3)/2*(TOL/(h0*sqrt(2))))*1/log(2));

% First we calculate the h-vector containing all the stepsizes for each
% level h(1) = h0, h(l) = h0*2^(-l+1), l>1
h = zeros(1, Lmax+1);
for l=1:Lmax+1
    h(l) = h0*2^(-l+1);
end

% Now we calculate the finest grid, for that we need the number of points
I = ceil(1./h);
I = [I 2/(h(end))];
x = linspace(0, 1, I(end)+1);
x_auxiliary = x(2:I(end));


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

% Check method
if strcmp(method, 'MC')
    % TODO Maybe add a shift, s.t. we can use Importance sampling
    Y = normrnd(0, 1, [cutoff M0]);
elseif strcmp(method, 'QMC')
    % TODO
    Y = normrnd(0, 1, [cutoff M0]);
elseif strcmp(method, 'RQMC')
    % TODO
    Y = normrnd(0, 1, [cutoff M0]);
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
    QoI = zeros(M0, 1);
    QoI_temp = zeros(M0, 1);
    for m=1:M0
        tic
        A = zeros(I(l)-1, I(l)-1);
        F = zeros(I(l)-1, 1);
        for i=1:I(l)-1
            F(i) = f(x(i));
            if i>1
                A(i, i-1) = - randomfield_a(indices(i), m)/h(l)^2;
            end
            A(i, i) = (randomfield_a(indices(i), m) + randomfield_a(indices(i+1), m))/h(l)^2;
            if i<I(l)-1
                A(i, i+1) = -randomfield_a(indices(i+1), m)/h(l)^2;
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
for l=1:Lmax
    if M_l(l) - M_l(l+1)<0
        error("Optimal sample size for each level is not decreasing.")
    end
end
M_l

% ------------------------------------------------------------------------
% Now we will calculate a MLMC estimate using the optimal number of samples
% we just calculated

% We can easily implement this procedure using the previous code with
% differnt sample sizes dependend on the level currently worked on

for l=1:Lmax+1
    % Check method and sample accordingly and corresponding to the current
    % level
    if strcmp(method, 'MC')
        % TODO Maybe add a shift, s.t. we can use Importance sampling
        Y = normrnd(0, 1, [cutoff M_l(l)]);
    elseif strcmp(method, 'QMC')
        % TODO
        Y = normrnd(0, 1, [cutoff M_l(l)]);
    elseif strcmp(method, 'RQMC')
        % TODO
        Y = normrnd(0, 1, [cutoff M_l(l)]);
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
    QoI = zeros(M_l(l), 1);
    for m=1:M_l(l)
        A = zeros(I(l)-1, I(l)-1);
        F = zeros(I(l)-1, 1);
        for i=1:I(l)-1
            F(i) = f(x(i));
            if i>1
                A(i, i-1) = - randomfield_a(indices(i), m)/h(l)^2;
            end
            A(i, i) = (randomfield_a(indices(i), m) + randomfield_a(indices(i+1), m))/h(l)^2;
            if i<I(l)-1
                A(i, i+1) = -randomfield_a(indices(i+1), m)/h(l)^2;
            end
        end
        u_h = A\F;

        if l == 1
            % QoI for level 0: sampling g0
            QoI(m) = h(l)*sum(u_h);
        else
            % QoI for level l>0: sampling g_l - g_{l-1}
            QoI(m) = h(l)*sum(u_h) - QoI_previous_level(m);
        end
    end

    % Calculate MC estimate, variance, avg. cost per sample for QoI
    MC_estimate(l) = 1/M_l(l)*sum(QoI);
    variance(l) = 1/(M_l(l)-1)*sum((QoI - MC_estimate(l)).^2);
    QoI_previous_level = QoI;
end

% Finally we can write the MLMC estimator
MLMC_estimate = sum(MC_estimate(1:Lmax));
MLMC_variance = sum((1./M_l(1:Lmax)).*variance(1:Lmax));

end