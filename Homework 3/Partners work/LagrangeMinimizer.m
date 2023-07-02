function [msol,Msol] = LagrangeMinimizer(V,delta,TOL,C)
    fun = @(M) M(8)*(C(1)*M(1)+C(2)*M(2)+C(3)*M(3)+C(4)*M(4)+C(5)*M(5)+C(6)*M(6)+C(7)*M(7));
    
    nonlcon = @(M) con(M,delta,TOL,V);
    
    A=-eye(8,8);
    b=zeros(8,1);
    Aeq=[];
    beq=[];
    lb=[];
    ub=[];

    x0=[1000 100 50 30 20 10 5 1];
    sol = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
    msol=sol(8);
    Msol=sol(1:7);
    
end

function [c,ceq] = con(M,delta,TOL,V)
    c = 1/M(8)*(V(1)/M(1)^(1+delta)+V(2)/M(2)^(1+delta)+V(3)/M(3)^(1+delta)+V(4)/M(4)^(1+delta)+V(5)/M(5)^(1+delta)+V(6)/M(6)^(1+delta)+V(7)/M(7)^(1+delta));
    ceq = [];
end