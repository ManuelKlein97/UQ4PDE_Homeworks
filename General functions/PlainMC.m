function [MC_estimate, MC_variance] = PlainMC(samples)
%{
---------------------------------------------------------------------------
Description:
     
---------------------------------------------------------------------------
Inputs:
      samples: samples vector (length(samples) x 1)

Outputs:
     MC_estimate:
     MC_variance:
---------------------------------------------------------------------------
%}
M = length(samples);
MC_estimate = 1/M*sum(samples);
MC_variance = 1/(M-1)*sum((samples - MC_estimate).^2);
end