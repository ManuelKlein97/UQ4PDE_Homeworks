function R = GenerateSamplesIS(mu, N, M)
    R = i4_sobol_generate(N,M,0);
    R = R + randn(N ,M) + mu;
    R = mod(R, 1);
    R = norminv(R);
end