function [solution] = Adaptive_MLQMC(I, M0, M_increment, TOL, nu, K, shift, breakvalue)
%{
---------------------------------------------------------------------------
Description:
Calculates Adaptive Monte Carlo  Estimator and Variance for P(Q(u)<K) 
and a given number of mesh points I and number of samples M0 and 
adaptively updates the number of samples until a relative error
of TOL is given. The QoI we are interested in will be determined 
by a random field with parameters nu (0.5 or infty). 
---------------------------------------------------------------------------
Inputs:
             I: Number of mesh points
            M0: Starting value for the number of samples in use
   M_increment: Determines how many more samples each iterations are used
           TOL: Tolerance which is supposed to be achieved by the procedure
            nu: random field parameter
             K: K-Value for which the QoI is supposed to be lower
         shift: Possible Importance sampling shift (optimization step needs to
                be done prior to this procedure)
    breakvalue: Determines after which iteration to stop
        method: QMC or RQMC

Outputs:
            estimate: MC estimate vector (1:M) estimates
            variance: Variance vector (1:M)
                   M: Final number of steps used
          iterations: Number of iterations needed for reaching TOL
      relative_error: Total relative error in each iteration (bias+stat)
          bias_error: Relative Bias error in each iteration
          stat_error: Relative Statistical error in each iteration
---------------------------------------------------------------------------
%}

end