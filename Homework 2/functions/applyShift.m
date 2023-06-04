function R = applyShift(R,S)
    %R = i4_sobol_generate(N,M,0);
    R = R + S;
    R = mod(R,1);
    R = norminv(R);
end