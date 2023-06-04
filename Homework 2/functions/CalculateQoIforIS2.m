function [V1] = CalculateQoIforIS2(rho,sigma,c_alpha, I, h,nu,cutoff,Y,K,EWsort,EVsort)
    F1 = assemble_F(I,h,@f);
    
    A = assemble_A(@a2_2,I,h,Y,cutoff,EWsort,EVsort);
    V1 = h * sum(A\F1);
end