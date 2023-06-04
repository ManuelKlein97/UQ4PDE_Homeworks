function [c,ceq] = nzcon(rho,sigma,c_alpha,I,h,nu,x,k,EWsort,EVsort)
    c = CalculateQoIforIS2(rho,sigma,c_alpha, I, h,nu,5,x,k,EWsort,EVsort)-k;
    ceq = [];
end