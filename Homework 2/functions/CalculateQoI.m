function Q = CalculateQoI( I, h,fourier,cutoff,y,z)
    F1 = assemble_F(I,h,@f);
    
    A = assemble_A(@a2_3,I,h,fourier,cutoff,y,z);
    Q = h * sum(A\F1);
end