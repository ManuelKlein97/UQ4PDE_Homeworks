function [solution] = exact_ref_solutions(name, N)
%{
---------------------------------------------------------------------------
Description:
Returns the exact reference solution of the integrals of the functions f
---------------------------------------------------------------------------
Parameters:
    name: name of the function handle
    N: input dimension of function handle f
---------------------------------------------------------------------------
%}
if name == "Oscillator"
    c = 9/N*ones([1 N]);
    w = 1/2;

    inner_prod = 1./(1i.*c).*(exp(1i*c)-1);
    solution = real(exp(1i*2*pi*w)*prod(inner_prod));
end

if name == "Product Peak"
    c = 7.25/N*ones([1 N]);
    w = 1/2*ones([1 N]);
    inner_prod = ones([1 N]);
    for l=1:N
        inner_prod(l) = c(l)*(atan(c(l)*(1-w(l))) + atan(c(l)*w(l)));
    end
    solution = prod(inner_prod);
end

if name == "Gaussian"
    inner_prod = ones([1 N]);
    w = 1/2;
    c = (7.03/N)*ones([1 N]);

    for l=1:N
        inner_prod(l) = sqrt(pi)/(2*c(l))*(erf(c(l)*(1-w)) + erf(c(l)*w));
    end
    solution = prod(inner_prod);
end

if name == "Continuous"
    w = 1/2;
    c = (2.04/N)*ones([1 N]);
    inner_prod = ones([1 N]);
    for l=1:N
        inner_prod(l) = 1/(c(l))*(2 - exp(-c(l)*w) - exp(-c(l)*(1-w)));
    end
    solution = prod(inner_prod);
end

if name == "Discontinuous"
    w = [pi/4 pi/5];
    c = (4.3/N)*ones([1 N]);
    solution = 1/prod(c)*(exp(c(1)*w(1)) - 1)*(exp(c(2)*w(2)) - 1)*prod(exp(c(3:N)) - 1);
end

end