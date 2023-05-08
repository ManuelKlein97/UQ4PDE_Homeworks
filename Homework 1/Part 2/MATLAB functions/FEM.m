function [solution] = FEM(I, randomfield)
%{
---------------------------------------------------------------------------
Description:
Function to return the approximate solution u_h of the given PDE
---------------------------------------------------------------------------
Parameters:
             I: Grid size for the uniform mesh to discret. problem
   randomfield: 1 == Piecewise constant coeff, 2 == Log-Normal
---------------------------------------------------------------------------
%}

% Check if input string for random field is correct
if ~(randomfield == 1 | randomfield == 2)
    error("Please hand over a value 'CCP'==1 or 'Log-Normal'==2 for the choice of the random field")
end


% Initialize uniform grid
x = linspace(0, 1, I+1);

% Set step size h
h = x(2) - x(1);

% Initialize matrix A (discr. of the random field a) but only the
% diagonalentries as own vectors (to later use tridiag to solve system)
A = zeros([2*I+1 2*I+1]);
A_up = zeros([2*I 1]);
A_diag = zeros([2*I+1 1]);
A_low = zeros([2*I 1]);

% Generate the array of intermediate values
x_half = linspace(0, 1, 2*I+1);

% Initialize vector F
F = zeros([2*I-1 1]);
for i=1:2*I-1
    F(i) = 4*pi^2*cos(2*pi*x_half(i));
end

% Generate random field sample for all the values of x_half
if (randomfield == 1)
    N = 10+1; % Set N=10 for the mesh to define the random field
    sigma = 0.5; % Set sigma=0.5 for now
    a = a_constant_coeff(x_half, N, sigma);
end
if (randomfield == 2)
    nu = 0.5; % Set nu=0.5 for now
    %l = 3;
    %N = I*2^(-l);
    %a = a_KL_expansion(x, nu, N);
    a = a_matern(x_half, nu);
end


for i=2:2*I
    A(i, i-1) = -a(i-1)/h^2;
    A_low(i-1) = -a(i-1)/h^2;
    A(i, i) = (a(i-1) + a(i+1))/h^2;
    A_diag(i) = (a(i-1) + a(i+1))/h^2;
    A(i, i+1) = -a(i+1)/h^2;
    A_up(i) = -a(i+1)/h^2;
end

A = A(2:2*I, :);
solution = A\F;

end