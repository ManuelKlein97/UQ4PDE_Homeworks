function sng = standardNormal(y)
    sum=0;
    for i=1:size(y,2)
        sum=sum+y(i)^2;
    end
    sng=1/sqrt(2*pi)*exp(-0.5*sum);
end