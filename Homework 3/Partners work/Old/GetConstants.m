function [c, v] = GetConstants(f,d,M0)
    N = 100;
    c_vec = zeros(1,N);
    for n = 1:N
        tic;
        Y=randn(d,1);
        f(Y);
        c_vec(n)=toc;
    end
    c = mean(c_vec);

    Y=randn(d,M0);
    temp=0;
    for i=1:M0
        temp=temp+f(Y(:,i));
    end
    temp=temp/M0;
    v=0;
    for i=1:M0
        v=v+(temp-f(Y(:,i)))^2;
    end
    v=v/(M0-1);
end


