function y = C(x,y,nu)
    if nu == 0.5
        y = 2 .* exp( - abs(x-y)/0.1 );
    end
    if nu == 1.5
        y = 2 .* (1 + sqrt(3)*abs(x-y)/0.1) .* exp( - sqrt(3)*abs(x-y)/0.1 );
    end
    if nu == 2.5
        y = 2 .* (1 + sqrt(5)*abs(x-y)/0.1 + sqrt(3)*abs(x-y).^2/0.1^2) .* exp( - sqrt(5)*abs(x-y)/0.1 );
    end
    if nu == 0 % for infinity
        y = 2 .* exp( - abs(x-y).^2/(2*0.1^2) );
    end
end