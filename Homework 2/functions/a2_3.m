function a = a2_3(x,I,fourier,cutoff,y,z)
    N = cutoff;
    if I < cutoff
        N=I;
    end
   
    kappa = 0;
    for n = 1:N
           kappa = kappa + fourier(n)*( y(n)*cos(n*pi*x/2) + z(n)*sin(n*pi*x/2) );
    end
    a = exp(kappa);
end