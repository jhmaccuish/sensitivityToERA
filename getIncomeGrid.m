function [Ygrid, Q, minInc, maxInc, AIMEgrid] = getIncomeGrid

% A function that returns:
% 1. an income grid
% 2. A Markovian transition matrix (Q) over income realisations
% 3. A vector of minimum incomes in each year
% 4. A vector of maximum incomes in each year

%% ------------------------------------------------------------------------ 
%  Declare the global variables that will be necessary in this function
global T Tretire sigma mu rho        % Structural economic paramters
global normBnd numPointsY                        % Numerical methods
global pension delta benefit numAIME %spouseInc
%% ------------------------------------------------------------------------ 

         sig_inc = sigma/((1-rho^2)^0.5);



%----------------------------------------%
%Now get a matrix, T * numPointsY that holds the grid for each income in
%each year. Do likewise with minimum and maximum income
%----------------------------------------%

%fc = table2array(readtable('1bias.txt','Delimiter','\r'));
fc=zeros(T,1);

AIMEgrid = NaN(T+1,numAIME);

[ly,Q] = tauchen(numPointsY,mu,rho,sigma,3);
Ygrid = exp(repmat(ly', T, 1)+repmat(polyval(delta,1:T)'-fc(1:T),1,numPointsY));
minInc = exp(repmat((-normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((-normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));          
maxInc = exp(repmat((normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));

pension = 107.45*52;
upper = cumsum(maxInc);
for t=1:T+1
    if t <=45
        AIMEgrid(t,:) = linspace(0,upper(t)/t,numAIME);
    else
        AIMEgrid(t,:) = AIMEgrid(t-1,:);
    end
end    

benefit = [57.90*52*ones(5,1);73.10*52*ones(length(minInc)-5,1)];%(73.10/107.45)*ones(length(minInc),1)*pension;%0.1*(minInc);

end