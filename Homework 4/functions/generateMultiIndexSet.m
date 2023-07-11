function [indxset] = generateMultiIndexSet(N, w, type)
%{
---------------------------------------------------------------------------
Description:
Function to return a multiindex set corresponding to a type used for N
dimensions and truncation parameter w
---------------------------------------------------------------------------
Parameters:
    N: dimension of space
    w: Truncation parameter for max polynomial degree
 type: Tensor product 'TP' / Hyperbolic Cross 'HC'
---------------------------------------------------------------------------
%}
if N == 1
    error('N is supposed to be greater than 1.')
end
% Tensor product
if strcmp(type, 'TP')    
    combinations = cell(1, N);
    values = 0:w;
    [combinations{:}] = ndgrid(values);
    indxset = reshape(cat(N+1, combinations{:}), [], N);
end

% Hyperbolic cross
if strcmp(type, 'HC')
    values = 0:w;
    % Generate all combinations with product less or equal to w
    combinations = combvec(values, values).';
    indxset = combinations(prod(combinations + 1, 2) <= w + 1, :);
end

% Total degree
if strcmp(type, 'TD')
    values = 0:w;
    % Generate all combinations with sum at most w
    combinations = combvec(values, values).';
    indxset = combinations(sum(combinations, 2) <= w, :);
end

end