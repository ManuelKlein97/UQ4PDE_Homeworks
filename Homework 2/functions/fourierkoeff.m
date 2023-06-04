function f = fourierkoeff(nu, cutoff)
f = zeros(1,cutoff);
    for n = 1:cutoff
        f(n) = sqrt(integral( @(x) C(x,0,nu).* cos(n*pi*x/2) ,0,2));
    end
end