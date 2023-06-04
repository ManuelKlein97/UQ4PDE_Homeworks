function y =  assemble_F(I,h,f)
    y = zeros(I-1,1);
    for i = 1:I-1
        y(i) = f( i*h);
    end
end