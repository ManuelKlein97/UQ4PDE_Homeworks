function [EWsort,EVsort] =  eigC(I,h,nu)
    sigma_a2 = sqrt(2);
    rho = 0.1;
    C_matrix = zeros(I,I);
        for m = 0:I-1
            for n = 0:I-1
                C_matrix(m+1,n+1) = C( ( m + 1/2)*h , ( n + 1/2)*h ,rho,sigma_a2,nu);
            end
        end
        [EV,EW] = eig(h*C_matrix);
    
        for m = 1:I
            EV(:,m) = EV(:,m) ./ (sqrt(h)*norm(EV(:,m) )); %norming the EV
        end
        EVsort = zeros(I,I);
        for m = 1:I
            EVsort(:,m) = EV(:,I-m+1);
        end
        EWsort = zeros(I,I);
        for m = 1:I
            EWsort(m,m) = EW(I-m+1,I-m+1);
        end
end