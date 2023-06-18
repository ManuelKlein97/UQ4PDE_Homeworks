function q = CalculateQoI(l,L,y,I,h,nu)%Diese Funktion bestimmt für einen festen Wert der ZV y einen Wert von Q
    sigma_a2 = sqrt(2);
    rho = 0.1;

    A = assemble_A_deterministic(@a2_2_deterministic,l,L,rho,sigma_a2,nu,I,h,y);
    F = assemble_F(l,L,I,h,@f);
    Q_h = A\F;
    q = h * sum(Q_h);
end
