function [c, v] = GetConstants(f, d, M0)
    % Returns time c that it takes to generate a sample and evaluate the
    % function f at the sample
    % v is the variance
    % we need both c and v for finding the optimal sample size in MLMC
    
    tic;
    Y=randn(d,1);
    f(Y);
    c=toc;

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


