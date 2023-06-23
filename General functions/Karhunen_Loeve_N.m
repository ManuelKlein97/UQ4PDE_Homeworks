function [randomfield] = Karhunen_Loeve_N(Y, EVal, EVec, N)
%{
---------------------------------------------------------------------------
Description:
    Returns a randomfield sample given randomsamples Y and Eigenvectors,
    Eigenvalues for a given model at a truncation level N.
---------------------------------------------------------------------------
Inputs:
        Y: Random variable samples (N0 x M) [M = Number of samples]
     EVal: Eigenvalues of a underlying model (N0 x 1)
     EVec: Corresponding Eigenvectors of the underlying model (N0 x (I-1))
        N: truncation scheme (N<=N0)
Outputs:
    randomfield: Randomfield sample using truncated N-term 
                 Karhunen Loeve Expansion, size (size_EVec x size(Y)(2)) or
                 
---------------------------------------------------------------------------
%}
[size_Evec, N1] = size(EVec);
I = size_Evec + 1;
[N2, M] = size(Y);

if N1<N || N2<N
    error("Truncation value N greater than size of samples.")
end

kappa = zeros(I-1, M);
randomfield = zeros(I-1, M);
for m=1:M
    for n=1:N
        kappa(:, m) = kappa(:, m) + sqrt(EVal(n))*Y(n, m)*EVec(:, n);
    end
    randomfield(:, m) = exp(kappa(:, m));
end

end