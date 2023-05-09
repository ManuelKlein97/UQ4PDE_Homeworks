function [solution] = crop(a)
%{
---------------------------------------------------------------------------
Description:
function to return a vector from a matrix a
---------------------------------------------------------------------------
Parameters:
         a: matrix with only zeros and one entry per row
---------------------------------------------------------------------------
%}

[n, m] = size(a);
boolarr = zeros([n 1]);
% Check if there are any overlapping entries in the rows
for k=1:n
    for l=1:m
        if boolarr(k) == 1
            if a(k, l) ~= 0
                a(k, l) = 0;
                %error("Double row entries.")
            end
        end
        if a(k, l) ~= 0
            boolarr(k) = 1;
        end

    end
end

solution = sum(a, 2);

end