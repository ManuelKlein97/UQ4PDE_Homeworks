function q = CalculateQoI(y,I,h,nu)%Diese Funktion bestimmt für einen festen Wert der ZV y einen Wert von Q
    sigma_a2 = sqrt(2);
    rho = 0.1;

    A = assemble_A_deterministic(@a2_2_deterministic,rho,sigma_a2,nu,I,h,y);
    F = assemble_F(I,h,@f);
    Q_h = A\F;
    q = h * sum(Q_h);
end


function y = assemble_A_deterministic(a,rho,sigma,nu,I,h,detVal)
    y = zeros(I-1,I-1);
    for i = 2:I-2
        x_right = (  (i-1)*h  +  i*h  )/2;
        x_left = ( (i-1)*h  + (i-2)*h )/2;
        [a_left,~,~] = a(x_left ,rho,sigma,nu,I,h,detVal);
        [a_right,~,~] = a(x_right ,rho,sigma,nu,I,h,detVal);
        y(i,i) = (a_left + a_right)/h^2;
        y(i,i+1) = -a_right/h^2;
        y(i,i-1) = -a_left/h^2;
    end
    [a_start,~,~] = a(h/2 ,rho,sigma,nu,I,h,detVal);
    [a_end,~,~] = a( ((I-1)*h + I*h)/2 ,rho,sigma,nu,I,h,detVal);
    [a_left,~,~] = a( ((I-1)*h + (I-2)*h)/2 ,rho,sigma,nu,I,h,detVal);
    [a_right,~,~] = a(3*h/2 ,rho,sigma,nu,I,h,detVal);
    y(1,1) = ( a_start + a_right )/h^2;
    y(I-1,I-1) = ( a_left  + a_end )/h^2;
    y(1,2) = -a_right/h^2 ;
    y(I-1,I-2) = -a_left/h^2;
end

function y =  assemble_F(I,h,f)
    y = zeros(I-1,1);
    for i = 1:I-1
        y(i) = f( i*h);
    end
end

function [a,EV,EW] = a2_2_deterministic(x,rho,sigma,nu,I,h,y)
    N = min(length(y),I);

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

    j = 0; 
    while j*h < x %the input x is always a multiple of h, so j will be that factor
        j=j+1;
    end
    
    kappa = 0;
    for n = 1:N
        kappa = kappa + sqrt( EW(n,n) )*y(n) * EV(j,n);
    end
    a = exp(kappa);
end

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

function y = f(x)
    y = 4*pi^2*cos(2*pi*x);
end

function y = C(x,y,rho,sigma,nu)
    if nu == 0.5
        y = sigma^2 .* exp( - abs(x-y)/rho );
    end
    if nu == 1.5
        y = sigma^2 .* (1 + sqrt(3)*abs(x-y)/rho) .* exp( - sqrt(3)*abs(x-y)/rho );
    end
    if nu == 2.5
        y = sigma^2 .* (1 + sqrt(5)*abs(x-y)/rho + sqrt(3)*abs(x-y).^2/rho^2) .* exp( - sqrt(5)*abs(x-y)/rho );
    end
    if nu == 0 % for infinity
        y = sigma^2 .* exp( - abs(x-y).^2/(2*rho^2) );
    end
end