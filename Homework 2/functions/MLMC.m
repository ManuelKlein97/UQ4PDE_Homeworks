function [MLMCres, MLMCvar, M_1, M_2] = MLMC(TOL, h_MC, h, Var_1, C_1, Var_2, C_2, nu,cutoff)
    for i = 1:length(h)
        if h(i) == h_MC
            index = i;
        end
    end
    
    M_1 = round((TOL/3)^(-2) * sqrt( Var_1(index)/mean(C_1(:,index)) ) * ( sqrt(Var_1(index)*mean(C_1(:,index))) + sqrt(Var_1(index)*mean(C_1(:,index))) ) );
    M_2 = round((TOL/3)^(-2) * sqrt( Var_2(index)/mean(C_2(:,index)) ) * ( sqrt(Var_1(index)*mean(C_1(:,index))) + sqrt(Var_2(index)*mean(C_2(:,index))) ) );

    h_use = h(index);
    I_use = 1/h_use;

    MCRes=zeros(2,1);
    y2 = GenerateSamples(cutoff,M_2);
    z2 = GenerateSamples(cutoff,M_2);

    y1 = GenerateSamples(cutoff,M_1);
    z1 = GenerateSamples(cutoff,M_1);

    V1 = zeros(1,M_2);
    V2 = zeros(1,M_2);
    V3 = zeros(1,M_1);
 
    fourier = fourierkoeff(nu,cutoff);
    for m2=1:M_2
        V1(m2) = CalculateQoI( I_use, h_use,fourier,cutoff,y2(:,m2),z2(:,m2));
        V2(m2) = CalculateQoI( 2*I_use, h_use/2,fourier,cutoff,y2(:,m2),z2(:,m2));

        MCRes(2,1) = MCRes(2,1) + V1(m2) - V2(m2);
    end

    for m1 = 1:M_1
        V3(m1) = CalculateQoI( I_use, h_use,fourier,cutoff,y1(:,m1),z1(:,m1));

        MCRes(1,1) = MCRes(1,1) + V3(m1);
    end
    MCRes(2,1) = MCRes(2,1)/M_2;
    MCRes(1,1) = MCRes(1,1)/M_1;
    Var=0;
    Var_diff = 0;
    for m2=1:M_2
        Var_diff = Var_diff + ( V1(m2) - V2(m2) -MCRes(2,1))^2;
    end
    for m1 = 1:M_1
        Var = Var + ( V3(m1) -MCRes(1,1) )^2;
    end
    MLMCres = MCRes(2,1) + MCRes(1,1);
    Var1 = Var/(M_1-1);
    Var2 = Var_diff/(M_2-1);    
    MLMCvar = Var1/M_1 + Var2/M_2;
end