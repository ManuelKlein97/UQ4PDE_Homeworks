function [solution] = QoI_parameter(I, randomfield, RF_par, N)
%{
---------------------------------------------------------------------------
Description:
Function to return the value of the QoI
---------------------------------------------------------------------------
Parameters:
             I: Grid size for the uniform mesh to discret. problem
   randomfield: 1 == Piecewise constant coeff, 2 == Log-Normal (KL exp.)
        RF_par: Random field parameter either nu or sigma
             N: 
             l:
---------------------------------------------------------------------------
%}

% Check if input string for random field is correct
if ~(randomfield == 1 | randomfield == 2)
    error("Please hand over a value 'CCP'==1 or 'Log-Normal'==2 for the choice of the random field")
end

if ~exist('RF_par', 'var')
    RF_par = 0.5;
end

if randomfield == 1
    if ~(RF_par < 1/sqrt(3))
        error("RF_par has to satisfy 1-RF_par*sqrt(3)>0 to ensure coercivity of the CCP model.")
    end
end

if randomfield == 2
    if ~(RF_par == 0.5 | RF_par == 1.5 | RF_par == 2.5 | RF_par == 'infinity')
        error("RF_par has to be either value 0.5, 1.5, 2.5 or string 'infinity' for the Log-Normal model.")
    end
end


% Initialize uniform grid
x = linspace(0, 1, I+1);

% Set step size h
h = x(2) - x(1);

% Initialize the x_{i+1/2} point array for evaluating RF for A construction
x_half = zeros([I+1 1]);
for k=1:I
    x_half(k) = (x(k+1) + x(k))/2; 
end


% Initialize vector F (rhs of PDE)
F = zeros([I+1 1]);
for i=1:I+1
    F(i) = 4*pi^2*cos(2*pi*x_half(i));
end

if (randomfield == 1)
    a = a_constant_coeff(x_half, N, RF_par);
end
if (randomfield == 2)
    a = a_KL_expansion(x_half, RF_par, N);
end

% Preallocate memory for tridiagonal matrix A
A = zeros([I+1 I+1]);
% Also for the diagonals only (to deploy efficient method for solving
% tridiagonal linear system)
A_up = zeros([2*I 1]);
A_diag = zeros([2*I+1 1]);
A_low = zeros([2*I 1]);

for i=2:I
   A(i, i-1) = -a(i-1)/(h^2);
   A(i, i) = (a(i-1) + a(i+1))/(h^2);
   A(i, i+1) = -a(i+1)/(h^2);
end

h = 1/I;
u_h = A(2:I, :)\F(2:I);
solution = h*sum(u_h(1:end-1));

end