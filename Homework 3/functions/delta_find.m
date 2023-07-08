function [delta] = delta_find(m, M)
V = [17.62 9.71 0.92 0.35 0.11 0.05 0.006];
C = [3.1571e-5 8.4513e-5 1.4618e-4 2.897e-4 6.7645e-4 1.3e-3 2.5e-3];
L = 7;

m1 = length(m);
m2 = length(M);
delta = zeros([m1 m2]);

for i=1:m1
    for j=1:m2
        delta(i, j) = 1/m(i)*sum(V./M(j));
    end
end


end