function [ value ] = objectivefunc(A1, A0, Y,L,ixP) % AIME

%%
%-------------------------------------------------------------------------------%
% This function returns the following quantity:
% - (u(c) +  b V( A1))
% where c is calculated from today's assets and tomorrow's assets
%-------------------------------------------------------------------------------%

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global beta r interpMethod            % structural model parameters
global Agrid1 EV1 AIME1grid                 % tomorrow's asset grid and tomorrow's expected value function

mortal = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.011071 0.011907 0.012807...
          0.013676 0.01475 0.015818 0.016846 0.018079 0.019343 0.020659 0.0218 0.023505 0.025202 0.02696...
          0.028831 0.031017 0.033496 0.036024 0.038912 0.042054 0.045689 0.049653 0.054036 0.05886 0.064093...
          0.069636 0.07533 0.081069 0.086912 0.093067 0.099807 0.107483 0.116125 0.125196 0.134361 0.143881...
          0.1542 0.165675 0.17842 0.192363 1];

%% ------------------------------------------------------------------------ 
%Get tomorrow's consumption (cons), the value of left over assets (VA1) and
%total value (u(c) + b * VA1
cons = A0  + Y - (A1)/(1+r);
VA1 = interp1(Agrid1,EV1,A1, interpMethod, 'extrap');
%VA1 = interp2(Agrid1, AIME1grid , EV1, A1, AIME);
%interp2(EV1, A1/Agrid1(20), AIME/AIME1grid(10));
value = utility(cons,L) + beta * (1- mortal(ixP))* VA1;

%% ------------------------------------------------------------------------ 
%The optimisation routine that we will use searches for the minimum of the
%function. We want the maximum. So we multiply out function here by -1 so
%that the optimiser will fill the minimum of the negative of our function,
%i.e. the maximum of our functino 

value = - value;

end