function [res] = MLMC(f,M,h0,hinc,L,d)
    %{
    ---------------------------------------------------------------------------
    Description:
    
    ---------------------------------------------------------------------------
    Inputs:
                 f: function
                 M:
                h0:
              hinc:
                 L:
                 d:
    
    Outputs:
               res:
    ---------------------------------------------------------------------------
    %}

    h=zeros(L,1);
    h(1)=h0;
    I=zeros(L,1);
    I(1)=1/h0;
    for i=2:L
        h(i)=h(i-1)/hinc;
        I(i)=I(i-1)*hinc;
    end
    temp=zeros(L,1);
    %M
    Y=randn(d,max(M));
    for i=1:M(1)
        temp(1)=temp(1)+f(h(1),Y(:,i));
    end
    temp(1)=temp(1)/M(1);

    for j=2:L
        %Y=randn(d,M(j));
        for i=1:M(j)
            temp(j)=temp(j)+f(h(j),Y(:,i))-f(h(j-1),Y(:,i));
        end
        temp(j)=temp(j)/M(j);
    end
    res=sum(temp);
end

