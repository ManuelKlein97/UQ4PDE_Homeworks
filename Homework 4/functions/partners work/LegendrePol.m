function LegrendrePol = LegendrePol(maxdeg)
    %{
    ---------------------------------------------------------------------------
    Description:
    Returns a matrix whose (i)-th column correspond to the coefficients of the
    (i-1)-th legendre polynomial
    ---------------------------------------------------------------------------
    Parameters:
          maxdeg: maximum degree of legendre polynomials displayed
        
    ---------------------------------------------------------------------------
    %}

    LegrendrePol = zeros(maxdeg,maxdeg); % The columns are the coefficients of the Legendre Polynomials.
    % In the first row is the coefficient for the leading exponent, so that
    % polyval(LegrendrePol(:,2) ,3) evaluates the second pol at x = 3
    P = diag(ones(1,maxdeg-1),1);% In the recursive definition one polynomial is multiplied by x, this is done with P
    I = eye(maxdeg);
    LegrendrePol(:,1) = I(:,maxdeg);
    LegrendrePol(:,2) =  I(:,maxdeg-1);
    for n = 3:maxdeg
        LegrendrePol(:,n) = (2*(n-1)-1)/(n-1)*P*LegrendrePol(:,n-1)  - ((n-1)-1)/(n-1)* LegrendrePol(:,n-2);
    end
    for n = 1:maxdeg
        LegrendrePol(:,n) = sqrt(2*(n-1)+1)* LegrendrePol(:,n);%norming the polynomials
    end
end