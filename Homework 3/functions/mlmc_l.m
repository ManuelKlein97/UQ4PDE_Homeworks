function [sums, costs] = mlmc_l(l, N0, varargin)

% inputs:  l = level
%          N = number of samples
%          varargin = optional additional user variables
%
% output: sums(1) = sum(Y)
%         sums(2) = sum(Y.^2)
%         where Y are iid samples with expected value:
%         E[P_0]           on level 0
%         E[P_l - P_{l-1}] on level l>0
%         cost = cost of N samples

% Noah: CalculateQoI(x,2I,1/(2I),nu) - CalculateQoI(x,I,1/I,nu)
%{
Goal: Compute the level l MC estimate for our specific problem (rPDE)
For that we need:
    1) A fixed Stepsize for level l=0,...,Lmin,...,Lmax (Discretization
    parameter hl=h0*2^(-l)
    2) 
    3)
%}


sums = [1 0.01];
costs = 0.001;

end