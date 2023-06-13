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

end