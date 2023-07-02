function [m] = lagrangian_MCQMC(V, C, delta, TOL)
L = length(V);
K = length(delta);
for k=1:K
    temp1 = (V.*(1+delta(k))./C).^(1./(1+delta(k))).*C;
    a = sum(temp1);
    
    %mu = 1/((1/(2+delta)) - 1 + (1+delta)/(2+delta))
    %nu = mu/((2+delta)*(1+delta)) - 2*mu/(2+delta) - 2/(2+delta)
    
    temp2 = V./((V.*(1+delta(k)))./C).^((1+delta(k))./(2+delta(k)));
    b = sum(temp2);
    %c = (b/a)^(mu/(2+delta))*((V*(1+delta))./C).^(1/(2+delta))
    
    m(k) = (a/b)^(4-1./(1+delta(k))-(2+2.*delta(k))./(2+delta(k)));
end
    %M_l = c*m.^nu
end