function [A] = getA(randomfield)
%{
---------------------------------------------------------------------------
Description:
     Returns a single realization of the random matrix A of the random PDE
     given a single randomfield sample
---------------------------------------------------------------------------
Inputs:
      randomfield: randomfield vector (length(randomfield) x 1)

Outputs:
     A: realization of random tridiagonal matrix corresponding to our rPDE
        (length(randomfield)-1 x length(randomfield)-1)
---------------------------------------------------------------------------
%}

% Goal: we need to find our the size of A depending on the randomfields
% size

temp = length(randomfield);
I = temp;
h = 1/(temp); % ok


% Now we need to get the underlying grid for the discretization
x_grid = grid_level(h); % ok

% not sure if we really need x_grid, bc we just index over RF with i=1:I


% Initialize Matrix A as sparse matrix
S = zeros(I, I);
A = sparse(S);

for i=1:I
    if i == 1
    % Case if h/2 has to be used
        A(i, i) = (randomfield(i) + randomfield(i+1))/((h/2)^2);
        A(i, i+1) = -randomfield(i+1)/((h/2)^2);
    elseif i == I
        A(i, i-1) = -randomfield(i)/((h/2)^2);
        A(i, i) = (randomfield(i-1) + randomfield(i))/((h/2)^2);
    else
        A(i, i-1) = -randomfield(i)/(h^2);
        A(i, i) = (randomfield(i) + randomfield(i+1))/(h^2);
        A(i, i+1) = -randomfield(i+1)/(h^2);
    end
end
end