function A = assemble_A(a,I,h,fourier,cutoff,y,z)
    A = zeros(I-1,I-1);
    for i = 2:I-2
        x_right = (  (i-1)*h  +  i*h  )/2;
        x_left = ( (i-1)*h  + (i-2)*h )/2;
        a_left = a(x_left ,I,fourier,cutoff,y,z);
        a_right = a(x_right ,I,fourier,cutoff,y,z);
        A(i,i) = (a_left + a_right)/h^2;
        A(i,i+1) = -a_right/h^2;
        A(i,i-1) = -a_left/h^2;
    end
    a_start = a(h/2 ,I,fourier,cutoff,y,z);
    a_end = a( ((I-1)*h + I*h)/2 ,I,fourier,cutoff,y,z);
    a_left = a( ((I-1)*h + (I-2)*h)/2 ,I,fourier,cutoff,y,z);
    a_right = a(3*h/2 ,I,fourier,cutoff,y,z);
    A(1,1) = ( a_start + a_right )/h^2;
    A(I-1,I-1) = ( a_left  + a_end )/h^2;
    A(1,2) = -a_right/h^2 ;
    A(I-1,I-2) = -a_left/h^2;
end