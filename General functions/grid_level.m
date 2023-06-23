function [grid, avg_grid] = grid_level(h)
%{
---------------------------------------------------------------------------
Description:
    Returns the grid for evaluating each level on interval [0, 1]
---------------------------------------------------------------------------
Inputs:
          h: stepsize

Outputs:
     grid: "Uniformly" spaced grid with I+1 points
 avg_grid: Averaged-out grid for evaluating Cov & RF
---------------------------------------------------------------------------
%}
I = 1/h;
start = h/2;
grid = [start];
for i=1:I-1
    grid = [grid start+i*h];
end
temp_grid = [grid 1];
temp_grid = [0 temp_grid];
avg_grid = zeros(I+1, 1);
for j=1:I+1
    avg_grid(j) = (temp_grid(j) + temp_grid(j+1))/2;
end

end