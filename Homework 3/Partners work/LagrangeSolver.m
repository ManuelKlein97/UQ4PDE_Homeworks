function [msol,Msol,lambdasol] = LagrangeSolver(V,delta,TOL)
%INPUT:
%V: Vector containing the variance at the different levels V_l
%delta: QMC constant
%TOL: wanted Tolerance
sol = sym('x', [9 1]);
options=optimoptions(@fsolve,'MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
optfun=@(x)f(x,V,delta,TOL);
assume(sol>=0)
sol=fsolve(optfun,[1,50,100000,500,250,125,100,50,10],options);
sol=real(sol);
lambdasol=sol(1);
msol=sol(2)^2;
Msol=zeros(7,1);
Msol(1)=sol(3)^2;
Msol(2)=sol(4)^2;
Msol(3)=sol(5)^2;
Msol(4)=sol(6)^2;
Msol(5)=sol(7)^2;
Msol(6)=sol(8)^2;
Msol(7)=sol(9)^2;
end

function y = f(x,V,delta,TOL)
    lambda=x(1);
    m=x(2)^2;
    M0=x(3)^2;
    M1=x(4)^2;
    M2=x(5)^2;
    M3=x(6)^2;
    M4=x(7)^2;
    M5=x(8)^2;
    M6=x(9)^2;
    y(1)=1/m*(V(1)/M0^(1+delta)+V(2)/M1^(1+delta)+V(3)/M2^(1+delta)+V(4)/M3^(1+delta)+V(5)/M4^(1+delta)+V(6)/M5^(1+delta)+V(7)/M6^(1+delta))-TOL^2;
    y(2)=M0+2*M1+4*M2+8*M3+16*M4+32*M5+64*M6-lambda/m^2*(V(1)/M0^(1+delta)+V(2)/M1^(1+delta)+V(3)/M2^(1+delta)+V(4)/M3^(1+delta)+V(5)/M4^(1+delta)+V(6)/M5^(1+delta)+V(7)/M6^(1+delta));
    y(3)=m-lambda*V(1)*(1+delta)/m*M0^(-2-delta);
    y(4)=2*m-lambda*V(2)*(1+delta)/m*M1^(-2-delta);
    y(5)=4*m-lambda*V(3)*(1+delta)/m*M2^(-2-delta);
    y(6)=8*m-lambda*V(4)*(1+delta)/m*M3^(-2-delta);
    y(7)=16*m-lambda*V(5)*(1+delta)/m*M4^(-2-delta);
    y(8)=32*m-lambda*V(6)*(1+delta)/m*M5^(-2-delta);
    y(9)=64*m-lambda*V(7)*(1+delta)/m*M6^(-2-delta);
end

