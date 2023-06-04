function [stepsize,samples,MCRes,StatisticalEr,BiasEr,RE,Variance_level,rt] = AdaptiveMC(M0, p,c_alpha, TOL, I0, h0,nu,cutoff)
    h=h0;
    rt=1;
    I=I0;
    MCRes=zeros(2,30);
    M=M0;
    QDiff=zeros(1,30);

    y = GenerateSamples(cutoff,M);
    z = GenerateSamples(cutoff,M);
    Bias = 0;

    fourier = fourierkoeff(nu,cutoff);
    V1 = zeros(1,M);
    V2 = zeros(1,M);
    for m=1:M
        V1(m) = CalculateQoI( I, h,fourier,cutoff,y(:,m),z(:,m));
        V2(m) = CalculateQoI( 2*I, h/2,fourier,cutoff,y(:,m),z(:,m));
 
        MCRes(1,rt)=MCRes(1,rt)+V1(m);
        MCRes(2,rt)=MCRes(2,rt)+V2(m);

        Bias=Bias+V1(m)-V2(m);
    end
    QDiff(rt)=abs(Bias)/M;
    Bias=QDiff(rt)/(1-0.5^p);
    MCRes(:,rt)=MCRes(:,rt)/M;
    Var=0;
    for m=1:M
        Var=Var+(V1(m)-MCRes(1,rt))^2;
    end
    Var=Var/(M-1);
    Bias = abs(Bias/MCRes(2,rt));
    stat = sqrt(Var)*c_alpha/(sqrt(M)*abs(MCRes(2,rt)));
    RE=[Bias+stat];

    if Bias > 0.5*TOL
		h=h/2;
	    I=2*I;
    end
		
    if stat > 0.5*TOL
	    M=2*M;
    end

    Variance_level = [Var];
	BiasEr = [Bias];
	StatisticalEr = [stat];
    rt=rt+1;
    while ((Bias > 0.5*TOL) || (stat > 0.5*TOL))
        if rt == 30
            break
        end
        Bias=0;
        y = GenerateSamples(cutoff,M);
        z = GenerateSamples(cutoff,M);
        V1 = zeros(1,M);
        V2 = zeros(1,M);
        for m=1:M
            V1(m) = CalculateQoI( I, h,fourier,cutoff,y(:,m),z(:,m));
            V2(m) = CalculateQoI( 2*I, h/2,fourier,cutoff,y(:,m),z(:,m));
     
            MCRes(1,rt)=MCRes(1,rt)+V1(m);
            MCRes(2,rt)=MCRes(2,rt)+V2(m);
    
            Bias=Bias+V1(m)-V2(m);
        end
        QDiff(rt)=Bias/M;
        Bias=QDiff(rt)/(1-0.5^p);
        MCRes(:,rt)=MCRes(:,rt)/M;
        Var=0;
        for m=1:M
            Var=Var+(V1(m)-MCRes(1,rt))^2;
        end
        Var=Var/(M-1);

        Bias = abs(Bias/MCRes(2,rt));
        stat = sqrt(Var)*c_alpha/(sqrt(M)*abs(MCRes(2,rt)));
        RE=[RE Bias+stat];
        Variance_level = [Variance_level Var];
		if Bias > 0.5*TOL
			BiasEr = [BiasEr Bias];
		    h=h/2;
			I=2*I;
        end
		if stat > 0.5*TOL
			StatisticalEr = [StatisticalEr stat];
			M=2*M;
        end
		rt=rt+1;
    end
    BiasEr = [BiasEr Bias];
    StatisticalEr = [StatisticalEr stat];
    Variance_level = [Variance_level Var];
    samples=zeros(length(StatisticalEr),1);
    samples(1)=M0;
    for m=2:length(StatisticalEr)
        samples(m)=samples(m-1)*2;
    end
	stepsize=zeros(length(BiasEr),1);
    stepsize(1)=h0;
    for m=2:length(BiasEr)
        stepsize(m)=stepsize(m-1)/2;
    end
end