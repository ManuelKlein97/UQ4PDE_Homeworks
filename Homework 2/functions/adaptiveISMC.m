function [estimate, variance, M, iterations, relative_error, bias_error, stat_error] = adaptiveISMC(I, M0, M_increment, TOL, nu, K, shift, breakvalue, method)
%{
---------------------------------------------------------------------------
Description:
Calculates Adaptive Monte Carlo Estimator and Variance for P(Q(u)<K) 
and a given number of mesh points I and number of samples M0 and 
adaptively updates the number of samples until a relative error
of TOL is given. The QoI we are interested in will be determined 
by a random field with parameters nu (0.5 or infty). 
---------------------------------------------------------------------------
Inputs:
             I: Number of mesh points
            M0: Starting value for the number of samples in use
   M_increment: Determines how many more samples each iterations are used
           TOL: Tolerance which is supposed to be achieved by the procedure
            nu: random field parameter
             K: K-Value for which the QoI is supposed to be lower
         shift: Possible Importance sampling shift (optimization step needs to
                be done prior to this procedure)
    breakvalue: Determines after which iteration to stop
        method: QMC or RQMC

Outputs:
            estimate: MC estimate vector (1:M) estimates
            variance: Variance vector (1:M)
                   M: Final number of steps used
          iterations: Number of iterations needed for reaching TOL
      relative_error: Total relative error in each iteration (bias+stat)
          bias_error: Relative Bias error in each iteration
          stat_error: Relative Statistical error in each iteration
---------------------------------------------------------------------------
%}

% Statistical part initialization
counter = 1;
fprintf('Iteration %s\n', num2str(counter));
M = M0;
h = 1/I;
c_alpha = 1.96;
QoI = zeros(1, M);
QoI_less_K = zeros(1, M);
QoI_double = zeros(1, M);
QoI_less_K_double = zeros(1, M);
estimate = zeros(1, M);
estimate_double = zeros(1, M);
variance = zeros(1, M);

% FEM initialization
F = zeros(I-1, 1);
F_double = zeros(2*I-1, 1);

% N-truncations of the fourier representation for nu=0.5 and nu=inf
% Cutoff should only be greater than 5 when there is no IS shift (i.e.
% shift == zeros)

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

if strcmp(nu, 'infinity')
    nu = 0;
end

%cutoff = length(shift);

% Generate normal samples (2 times)
sobolpoints1 = i4_sobol_generate(cutoff, M, 0);
rng(1); % Seed for reproducability
if strcmp(method, 'QMC')
    y = norminv(sobolpoints1, shift', 1);
elseif strcmp(method, 'RQMC')
    randomized_sobolpoints1 = sobolpoints1 + randn(cutoff, M) + shift';
    uniformed_rand_sp1 = mod(randomized_sobolpoints1, 1);
    y = norminv(uniformed_rand_sp1);
end

sobolpoints2 = i4_sobol_generate(cutoff, M, 0);
rng(2); % Seed for reproducability
if strcmp(method, 'QMC')
    z = norminv(sobolpoints2, shift', 1);
elseif strcmp(method, 'RQMC')
    randomized_sobolpoints2 = sobolpoints2 + randn(cutoff, M) + shift';
    uniformed_rand_sp2 = mod(randomized_sobolpoints2, 1);
    z = norminv(uniformed_rand_sp2);
end

% Now we can represent the random field a using the Fourier series
fourier = fourierkoeff(nu, cutoff);

% Solve problem for M once (after that comes a while loop)
% Solve FEM linear system of eq.

% Construction of F for normal mesh size
for i=1:I-1
    F(i) = f(i*h);
end

% Construction of F for double the mesh size
for i=1:2*I-1
    F_double(i) = f(i*h/2);
end

for m=1:M
    % Construction of A for normal mesh size
    A = assemble_A(@a2_3, I, h, fourier, cutoff, y(:, m), z(:, m));
    A_double = assemble_A(@a2_3, 2*I, h/2, fourier, cutoff, y(:, m), z(:, m));
    
    uh = A\F;
    uh_double = A_double\F_double;

    QoI(m) = h * sum(uh);
    

    if (m == M-1)
        QoI(m) = K - 0.1;
    end

    QoI_less_K(m) = QoI(m) < K;
    QoI_double(m) = h * sum(uh_double);
    QoI_less_K_double(m) = QoI_double(m) < K;
    estimate(m) = 1/m*sum(QoI_less_K(1:m));
    estimate_double(m) = 1/m*sum(QoI_less_K_double(1:m));
end

for m=1:M
    if m>1
        variance(m) = 1/(m-1)*sum((QoI_less_K(1:m) - estimate(end)).^2);
    end
end

% Check Relative Error (Statistical + Bias error) after M0 samples
%ref_sol = 0.0014;
stat_error = [sqrt(variance(M))*c_alpha/(sqrt(M)*abs(estimate_double(M)))];
if stat_error(counter) == 0
    stat_error(counter) = 1;
end

% TODO for bias estimation we need to calculate the QoI for double the mesh
% size 2*I
bias_error = [4/3*abs(1/M*(estimate(M) - estimate_double(M)))/abs(estimate_double(M))];
relative_error = [stat_error + bias_error];

% Now for the adaptive procedure
while relative_error(counter)>TOL

    % Set new M
    M = M + M_increment;
    % Generate normal samples (2 times)
    sobolpoints1 = i4_sobol_generate(cutoff, M, 0);
    rng(1); % Seed for reproducability
    if strcmp(method, 'QMC')
        y = norminv(sobolpoints1, shift', 1);
    elseif strcmp(method, 'RQMC')
        randomized_sobolpoints1 = sobolpoints1 + randn(cutoff, M) + shift';
        uniformed_rand_sp1 = mod(randomized_sobolpoints1, 1);
        y = norminv(uniformed_rand_sp1);
    end
    
    sobolpoints2 = i4_sobol_generate(cutoff, M, 0);
    rng(2); % Seed for reproducability
    if strcmp(method, 'QMC')
        z = norminv(sobolpoints2, shift', 1);
    elseif strcmp(method, 'RQMC')
        randomized_sobolpoints2 = sobolpoints2 + randn(cutoff, M) + shift';
        uniformed_rand_sp2 = mod(randomized_sobolpoints2, 1);
        z = norminv(uniformed_rand_sp2);    
    end
    
    % Now we can represent the random field a using the Fourier series
    fourier = fourierkoeff(nu, cutoff);
    
    
    % Solve problem for M once (after that comes a while loop)
    % Solve FEM linear system of eq.
    
    % Construction of F for normal mesh size
    for i=1:I-1
        F(i) = f(i*h);
    end
    
    % Construction of F for double the mesh size
    for i=1:2*I-1
        F_double(i) = f(i*h/2);
    end
    
    for m=1:M
        % Construction of A for normal mesh size
        A = assemble_A(@a2_3, I, h, fourier, cutoff, y(:, m), z(:, m));
        A_double = assemble_A(@a2_3, 2*I, h/2, fourier, cutoff, y(:, m), z(:, m));
        
        uh = A\F;
        uh_double = A_double\F_double;
    
        QoI(m) = h * sum(uh);
        
        
        if (m == M-1)
            QoI(m) = K - 0.1;
        end
        
        
        QoI_double(m) = h * sum(uh_double);
        
        QoI_less_K(m) = (QoI(m) < K);
        QoI_less_K_double(m) = (QoI_double(m) < K);

        estimate(m) = 1/m*sum(QoI_less_K(1:m));
        estimate_double(m) = 1/m*sum(QoI_less_K_double(1:m));
        
    end

    for m=1:M
        if m>1
            variance(m) = 1/(m-1)*sum((QoI_less_K(1:m) - estimate(end)).^2);
        end
    end


    stat_error = [stat_error sqrt(variance(end))*c_alpha/(sqrt(M)*abs(estimate_double(end)))];
    if stat_error(counter) == 0
        stat_error(counter) = 1;
    end

    bias_error = [bias_error 4/3*abs(1/M*(estimate(end) - estimate_double(end)))/abs(estimate_double(end))];
    relative_error = [relative_error stat_error(counter)+bias_error(counter)];
    
    fprintf("\nRel. err: %s", relative_error(counter))

    counter = counter + 1;
    fprintf('\nIteration %s\n', num2str(counter));
    if counter == breakvalue
        fprintf('\nForced stopping of procedure after %s iterations', num2str(breakvalue));
        break
    end
end
iterations = counter;
fprintf('\nRelative error: %s', num2str(relative_error(end)))
fprintf('\nStatistical error: %s', num2str(stat_error(end)))
fprintf('\nBias error: %s', num2str(bias_error(end)))

end