function [func_eval] = LegendreSeries(x, coeff)
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

[~, w] = size(coeff);
[~, n] = size(x);
func_eval = zeros(1, n);
for p=1:w
    func_eval = func_eval + coeff(p)*sqrt(2*(p-1)+1)*legendreP(p-1, 2*x-1);
end

end