function [sums,cost] = mlmc_l(l,N,varargin)
    if l>0
       f=@(x)(CalculateQoI(x,2^l*4,1/(2^l*4),0.5)-CalculateQoI(x,4*2^(l-1),1/(4*2^(l-1)),0.5));
    else 
       f=@(x)CalculateQoI(x,4,1/4,0.5);
    end
    tic;
    Y=randn(78,N);
    sums=zeros(2,1);
    for i=1:N
       sums(1)=sums(1)+f(Y(:,i));
    end
    cost=toc;
    for i=1:N
       sums(2)=sums(2)+f(Y(:,i))^2;
    end
    
end

