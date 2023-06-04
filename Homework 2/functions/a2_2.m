function a = a2_2(x,I,Y,cutoff,EWsort,EVsort)
    h = 1/I;
    N = cutoff;
    if I < cutoff
        N=I;
    end
    
    j = 0; 
    while j*h < x %the input x is always a multiple of h, so j will be that factor
        j=j+1;
    end

    kappa = 0;
    for n = 1:N
           kappa = kappa + sqrt( EWsort(n,n) )*Y(n) * EVsort(j,n);
    end
    a = exp(kappa);

end