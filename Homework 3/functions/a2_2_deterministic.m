function [a,EV,EW] = a2_2_deterministic(x,rho,sigma,nu,I,h,y) % Bestimmt das random field
    N = min(length(y),I);%y ist der rv, also ist length(y) der cutoff
    
    %EW und EV bestimmen
    C_matrix = zeros(I,I);
    for i = 0:I-1
        for j = 0:I-1
            C_matrix(i+1,j+1) = C( ( i + 1/2)*h , ( j + 1/2)*h ,rho,sigma,nu);
        end
    end
    [EV,EW] = eig(h*C_matrix);
    for i = 1:I
        EV(:,i) = EV(:,i) ./(sqrt(h)* norm(EV(:,i) )); %norming the EV
    end
    %EWsort und EVsort sind so sortiert, dass der EV mit größter spectraler
    %kontribution vorne steht
    EVsort = zeros(I,I); 
    for m = 1:I
        EVsort(:,m) = EV(:,I-m+1);
    end
    EWsort = zeros(I,I);
    for m = 1:I
        EWsort(m,m) = EW(I-m+1,I-m+1);
    end
    j = 0; 
    while j*h < x %the input x is always a multiple of h, so j will be that factor
        j=j+1;
    end
    %Jetzt kann das eigentliche rf bestimmt werden.  Das j in EVsort(j,n)
    %entspricht der deterministischen x-Koordinate
    kappa = 0;
    for n = 1:N
        kappa = kappa + sqrt( EWsort(n,n) )*y(n) * EVsort(j,n);
    end
    a = exp(kappa);
end
