function x = tridiag_solve(a,b,c,f)
% Solves the tridiagonal linear system Ax=f using the Thomas algorithm
% where A is a tridiagonal matrix with main diagonal b and subdiagonals a
% and c.

n = length(f);
x = zeros(n,1);

% Forward elimination
for i = 2:n
    w = a(i)/b(i-1);
    b(i) = b(i) - w*c(i-1);
    f(i) = f(i) - w*f(i-1);
end

% Back substitution
x(n) = f(n)/b(n);
for i = n-1:-1:1
    x(i) = (f(i)-c(i)*x(i+1))/b(i);
end

end