function [eigenvalues, eigenvectors] = getEigen(nu, cutoff, I)
%{
---------------------------------------------------------------------------
Description:
    Returns the eigenvalues and eigenvectors for a given parameter nu for
    model 2 until a certain cutoff.
---------------------------------------------------------------------------
Inputs:
        nu: Model 2 randomfield parameter
    cutoff: Truncation scheme
         I: Number of discretization points

Outputs:
     eigenvalues: size (cutoff x 1)
    eigenvectors: size (cutoff x I+1)
---------------------------------------------------------------------------
%}

rho = 0.1;
sigma = sqrt(2);
h = 1./I;

[temp, x_auxiliary] = grid_level(h);
NoP = length(x_auxiliary);

if cutoff > NoP
    cutoff = NoP;
end

Cov_mat = zeros(NoP, NoP);
for i=1:NoP
    for j=1:NoP
        if nu == 0.5
            Cov_mat(i, j) = sigma^2 .* exp( - abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif nu == 1.5
            Cov_mat(i, j) = sigma^2 .* (1 + sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j))/rho) .* exp( - sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif nu == 2.5
            Cov_mat(i, j) = sigma^2 .* (1 + sqrt(5)*abs(x_auxiliary(i)-x_auxiliary(j))/rho + sqrt(3)*abs(x_auxiliary(i)-x_auxiliary(j)).^2/rho^2) .* exp( - sqrt(5)*abs(x_auxiliary(i)-x_auxiliary(j))/rho );
        elseif strcmp(nu, 'infinity')
            Cov_mat(i, j) = sigma^2 .* exp(-abs(x_auxiliary(i)-x_auxiliary(j)).^2/(2*rho^2));
        end
    end
end
Cov_mat = h.*Cov_mat;

[unsorted_EVec, unsorted_EVal] = eig(Cov_mat);
for i = 1:NoP
    % Norm the EV
    unsorted_EVec(:, i) = unsorted_EVec(:, i)./(sqrt(h)*norm(unsorted_EVec(:, i)));
end
% Sort the Eigenvalues/Eigenvectors
EVec = zeros(NoP, NoP); 
EVal = zeros(NoP, 1);
for m = 1:NoP
    EVec(:, m) = unsorted_EVec(:, NoP+1-m);
    EVal(m) = unsorted_EVal(NoP+1-m, NoP+1-m);
end
eigenvalues = EVal(1:cutoff);
eigenvectors = EVec(:, 1:cutoff);

end