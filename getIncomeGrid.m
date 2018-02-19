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
global pension delta benefit numAIME spouseInc
%% ------------------------------------------------------------------------ 


% %----------------------------------------%
% % Scenario where there is no uncertainty
% %----------------------------------------%
% if isUncertainty == 0
%     y = exp(mu);                 % income is set equal to the exp of the log mean
%     minInc = y;                 
%     maxInc = y;
%     Q = [1];                    % The transition matrix Q is simply a constant 1
%                                 % with prob 1 each period income is 1
% 
% %----------------------------------------%
% % Scenario where there is uncertainty - income draws are log normally distributed 
% %----------------------------------------%
%                                 
% elseif isUncertainty == 1
% 
% 
%         % First get the standard deviation of income (from sigma and rho)
         sig_inc = sigma/((1-rho^2)^0.5);
% 
%         % Split the entire normal distribution into numPointsY sections that
%         % are equiprobable. The output lNormDev gives the (numPointsY + 1) 
%         % points that bound the sections, the output ly gives the 
%         % (numPointsY) expected value in each section
%         [ lNormDev, ly ] = getNormDev(mu, sig_inc, normBnd, numPointsY );    
% 
%         %---------------------%
%         %Get transition matrix Q(i, j). The prob of income j in t+1
%         %conditional on income i in t
%         %---------------------%
%         Q = NaN(numPointsY, numPointsY);             %initialise the transition matrix
%         for i = 1:1:numPointsY
%             for j = 1:1:numPointsY
%                 hiDraw = lNormDev(j+1) - (1-rho)*mu - rho * ly(i); %highest innovation that will give us income j tomorrow
%                 loDraw = lNormDev(j)   - (1-rho)*mu - rho * ly(i); %lowest  innovation that will give us income j tomorrow
%                 Q(i,j) = stdnormcdf_manual(hiDraw/sigma) - ...
%                          stdnormcdf_manual(loDraw/sigma);
%             end %j
%             
%             %Each of the rows of Q should add up to 1. But
%             %due to the truncation of the normal distribution they will
%             %not. So we divide through by the sum of elements in the row
%             Q(i, :) = Q(i, :) ./ sum(Q(i, :));                
%         end %i
% 
%         %y = exp(ly);                      % Get y from log y
%         %minInc = exp(-normBnd * sig_inc); %Get the minimum income in each year
%         %minInc = (-normBnd * sig_inc); %Get the minimum income in each year
%         %maxInc = exp(normBnd * sig_inc);  %Get the maximum income in each year 
%         %maxInc = (normBnd * sig_inc);  %Get the maximum income in each year 
%         
% %         if (y(1) < 1e-4) || (y(numPointsY) > 1e5)
% %             warning('Combination of sigma and rho give a very high income variance. Numerical instability possible')
% %         end
%         
% end  % if isUncertainty == 0


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
% ly = [ly; zeros(numPointsY,1)];
% Q = [Q.*repmat(1-probU,1,numPointsY) ,diag(probU) ;
%      diag(probE),diag(1-probE) ];
% Ygrid = exp(repmat(ly', T, 1)+repmat(polyval(delta,1:T)'-fc(1:T),1,2*numPointsY));
% minInc = exp(repmat((-normBnd * sigma/((1-rho^2)^0.5))', T, 1)+repmat(polyval(delta,1:T),1,2*numPointsY));          
% maxInc = exp(repmat((normBnd * sigma/((1-rho^2)^0.5))', T, 1)+repmat(polyval(delta,1:T),1,2*numPointsY));
pension = 107.45*52;%0.17*mean(mean(Ygrid));
Ygrid(Tretire:end,:)=Ygrid(Tretire:end,:)+pension;
%lower = cumsum(minInc);
upper = cumsum(maxInc);
for t=1:T+1
    if t <=45
        AIMEgrid(t,:) = linspace(0,upper(t)/t,numAIME);
    else
        AIMEgrid(t,:) = AIMEgrid(t-1,:);
    end
end    

Ygrid(45:end,:)=Ygrid(45:end,:)+pension;
Ygrid = Ygrid + spouseInc; 
benefit = [57.90*52*ones(5,1);73.10*52*ones(length(minInc)-5,1)];%(73.10/107.45)*ones(length(minInc),1)*pension;%0.1*(minInc);
benefit(Tretire:end)=benefit(Tretire:end)+pension;%(142.70/73.10)*benefit(Tretire:end,:);
benefit(45:end)=benefit(45:end)+pension;
benefit = benefit + spouseInc;

end