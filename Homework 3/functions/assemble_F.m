function y =  assemble_F(l,L,I,h,f) %Stellt den Vektor für die Funktion f auf
    if l == L%in this case the grid is equidistant and nothing needs to change
        y = zeros(I-1,1);
        for i = 1:I-1
            y(i) = f( i*h);
        end
    end
    if l < L
        points = [];
        for i =0:I-1
            points = [points (2*i+1)*h/2];
        end
        %points
        y = zeros(I,1);
        for i = 1:I
            y(i) = f(points(i));
        end
    end
end