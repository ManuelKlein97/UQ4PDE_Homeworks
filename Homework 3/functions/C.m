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
