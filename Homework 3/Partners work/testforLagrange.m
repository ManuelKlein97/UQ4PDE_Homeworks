V=[17.629 9.7068 0.9202 7.5329 6.5521 0.0519 0.0060];
delta=0.4;
TOL=0.1;
%V=[17.629 9.7068 0.9202 7.5329];
%LagrangeMinimizer(V,delta,TOL)
[m,M,lambda]=LagrangeSolver(V,delta,TOL);