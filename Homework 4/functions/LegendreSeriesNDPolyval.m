function [func_eval] = LegendreSeriesNDPolyval(x, coeff, indxset)
%{
---------------------------------------------------------------------------
Description:
Calculates the function evaluations of a Legendre Series expansion given
the legendre coefficients and a linspace x (or just a single value)
---------------------------------------------------------------------------
Parameters:
        x: linspace ([0,1])
    coeff: Legendre Coefficients
---------------------------------------------------------------------------
%}

[w, ~] = size(coeff);
w1 = max(indxset);
LegPolynomCoeff = LegendrePol(w1(1)+1);
[N, n] = size(x);
func_eval = zeros(1, n);
for p=1:w
    temp = ones(1, n);
    for i=1:N
        temp = temp*sqrt(2*indxset(p, i)+1).*polyval(LegPolynomCoeff(:, indxset(p, i)+1), 2*x(i, :)-1);
    end
    func_eval = func_eval + coeff(p).*temp;
end

end