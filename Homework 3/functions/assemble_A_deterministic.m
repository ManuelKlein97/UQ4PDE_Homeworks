function y = assemble_A_deterministic(a,l,L,rho,sigma,nu,I,h,detVal) %Stellt die Matrix A auf. detVal ist der RV
    if l == L%in this case the grid is equidistant and nothing needs to change
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
        [a_start,~,~] = a( h/2 ,rho,sigma,nu,I,h,detVal);
        [a_end,~,~] = a( ((I-1)*h + 1)/2 ,rho,sigma,nu,I,h,detVal);
        [a_left,~,~] = a( ((I-1)*h + (I-2)*h)/2 ,rho,sigma,nu,I,h,detVal);
        [a_right,~,~] = a(3*h/2 ,rho,sigma,nu,I,h,detVal);
        y(1,1) = ( a_start + a_right )/h^2;
        y(I-1,I-1) = ( a_left  + a_end )/h^2;
        y(1,2) = -a_right/h^2 ;
        y(I-1,I-2) = -a_left/h^2;
    end
    if l < L
        points = [];
        for i =0:I-1
            points = [points (2*i+1)*h/2];
        end
        %points
        y = zeros(I,I);
        for i = 2:I-1
            x_right = (points(i) + points(i+1))/2;
            x_left = (points(i) + points(i-1))/2;
            [a_left,~,~] = a(x_left ,rho,sigma,nu,I,h,detVal);
            [a_right,~,~] = a(x_right ,rho,sigma,nu,I,h,detVal);
            y(i,i) = (a_left + a_right)/h^2;
            y(i,i+1) = -a_right/h^2;
            y(i,i-1) = -a_left/h^2;
        end
        [a_start,~,~] = a( points(1)/2 ,rho,sigma,nu,I,h,detVal);
        [a_end,~,~] = a( (points(length(points)) + 1)/2 ,rho,sigma,nu,I,h,detVal);
        [a_left,~,~] = a( (points(length(points)) + points(length(points)-1))/2 ,rho,sigma,nu,I,h,detVal);
        [a_right,~,~] = a( (points(1)+points(2))/2 ,rho,sigma,nu,I,h,detVal);
        y(1,1) = ( a_start + a_right )/h^2;
        y(I,I) = ( a_left  + a_end )/h^2;
        y(1,2) = -a_right/h^2 ;
        y(I,I-1) = -a_left/h^2;
    end
end