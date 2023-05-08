function [solution] = QoI(I, randomfield, RF_par)
%{
---------------------------------------------------------------------------
Description:
Function to return the value of the QoI
---------------------------------------------------------------------------
Parameters:
             I: Grid size for the uniform mesh to discret. problem
   randomfield: 1 == Piecewise constant coeff, 2 == Log-Normal (KL exp.)
        RF_par: Random field parameter either nu or sigma
---------------------------------------------------------------------------
%}

if ~exist('RF_par', 'var')
    RF_par = 0.5;
end

h = 1/I;
u_h = FEM_KL(I, randomfield, RF_par);
solution = h*sum(u_h(1:end-1));


end