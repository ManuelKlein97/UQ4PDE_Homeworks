function opt = findbeststart(k, j, I, h, nu)
    opt = zeros(5,1);
    ID = eye(5);
    temp = zeros(5, 1);
    for i=1:5
        ith_id = ID(:, i);
        [temp, value] = fminsearch(@(x) minimizer(opt + x*ith_id, k, j, I, h, nu), 1);
        opt(i) = temp(i);
    end
end