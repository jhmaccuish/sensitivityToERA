function [ policyA1, policyC, policyL, V, EV, EdU ] = solveValueFunction

%This function obtains the value function for each time period and
%the policy function (i.e. optimal next-period asset choice) for each time
%period. From there we can work the optimal consumption level.

%The approach taken is by backwards recursion. The optimisation each period
%is carried out using 'fminbnd'. This is an in-built optimiser in Matlab.
%The optimisation routine it uses is known as the 'golden search method'

global T r tol minCons
global numPointsA numPointsY Agrid Ygrid incTransitionMrx numPointsL
global Agrid1 EV1 numAIME AIMEgrid AIME1grid
global benefit Tretire pension spouseInc

%% ------------------------------------------------------------------------ 
% GENERATE MATRICES TO STORE NUMERICAL APPROXIMATIONS AND INITIATE AS NAN

m=1;
% Matrices to hold the policy, value and marginal utility functions 
 V        = NaN(T+1, numPointsA, m*numPointsY, numAIME);
 policyA1 = NaN(T,   numPointsA, m*numPointsY, numAIME);
 policyC  = NaN(T,   numPointsA, m*numPointsY, numAIME);  
 policyL  = NaN(T,   numPointsA, m*numPointsY, numAIME);  

%dU       = NaN(T,   numPointsA, m*numPointsY, numAIME);

%Matrices to hold expected value and marginal utility functions 
EV  = NaN(T+1, numPointsA, m*numPointsY, numAIME);
EdU = NaN(T,   numPointsA, m*numPointsY, numAIME);


%% ------------------------------------------------------------------------ 
%Set the terminal value function and expected value function to 0
EV(T + 1, :,:)  = 0;          % continuation value at T-1
V(T + 1,:,:) = 0; 
%% ------------------------------------------------------------------------ 
% SOLVE RECURSIVELY THE CONSUMER'S PROBLEM, STARTING AT TIME T-1 AND MOVING
% BACKWARDS TO ZERO, ONE PERIOD AT A TIME

for ixt=T:-1:1                               % Loop from time T-1 to 1
    Agrid1 = Agrid(ixt + 1, :);               % The grid on assets tomorrow
    AIME1grid = AIMEgrid(ixt + 1, :); 
    for ixAIME = 1:1:numAIME
%         if ixt > 65 && ixAIME >1
%             policyA1(ixt,:,:,ixAIME)=policyA1(ixt,:,:,1);
%             policyL(ixt,:,:,ixAIME)=policyL(ixt,:,:,1);
%             policyC(ixt, :, :,ixAIME) =policyC(ixt, :, :,1);
%             V(ixt, :, :,ixAIME) = V(ixt, :, :,1);             
%             EV(ixt, :, :,ixAIME)  = EV(ixt, :, :,1);
%             continue
%         end
        for ixA = 1:1:numPointsA                   % points on asset grid
            %Although doesn't recieve income still need to loop around
            %hypothetical income because participation doesn't effect
            %earning potential
            % STEP 1. solve problem at grid points in assets, income + labour choices
            % ---------------------------------------------------------        
            for ixY = 1:1:numPointsY               % points on income grid
                negV = inf;
                for ixL = 0:1:(numPointsL-1)           % points of income choice
                        % Value of income and information for optimisation
                        A    = Agrid(ixt, ixA);            % assets today
                        Y    = ixL*Ygrid(ixt, ixY)+ (~ixL)*benefit(ixt);%+(ixt>=Tretire)*pension; % income today
                        if ixL &&  ixt < 65 
                            AIME =  Ygrid(ixt, ixY)/ixt + AIME * (ixt-1)/ixt;
                        elseif ~ixL && ixt >= 65
                            AIME = AIMEgrid(ixt,ixAIME);
                            Y =   Y + dbPension(AIME);
                        end   
                        %AIME only depends on earned income so add spousal
                        %and pensions after calculating it
                        Y = Y + (ixt >= Tretire)*pension+spouseInc + (ixt >= 45)*pension;
                        lbA1 = Agrid(ixt + 1, 1);          % lower bound: assets tomorrow
                        ubA1 = (A + Y - minCons)*(1+r);    % upper bound: assets tomorrow
                        EV1  = squeeze(EV(ixt + 1,:, ixY,:));  % relevant section of EV matrix (in assets tomorrow)

                        % Compute solution 
                        if (ubA1 - lbA1 < minCons)         % if liquidity constrained
                            negVtemp = objectivefunc(lbA1, A, Y,ixL,ixt, AIME);
                            policyA1temp = lbA1;
                        else                               % if interior solution
                            [policyA1temp, negVtemp] = ...
                                fminbnd(@(A1) objectivefunc(A1, A, Y,ixL,ixt, AIME), lbA1, ubA1, optimset('TolX',tol));  
                        end % if (ubA1 - lbA1 < minCons)
                        if negVtemp < negV
                            negV = negVtemp;
                            policyA1(ixt,ixA,ixY,ixAIME)=policyA1temp;
                            policyL(ixt,ixA,ixY,ixAIME)=ixL;
                        end
                end
                % Store solution and its value
                policyC(ixt, ixA, ixY,ixAIME) = A + Y - policyA1(ixt, ixA, ixY,ixAIME)/(1+r);
                V(ixt, ixA, ixY,ixAIME)       = -negV; 
                %dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY),policyL(ixt,ixA,ixY));            
            end 

            % STEP 2. integrate out income today conditional on income
            % yesterday to get EV and EdU
            % --------------------------------------------------------
            realisedV(:,:) = V(ixt, ixA, :, ixAIME);
            %realiseddU(:,:) = dU(ixt, ixA, :);
            for ixY = 1:1:numPointsY
                EV(ixt, ixA, ixY, ixAIME)  = incTransitionMrx(ixY,:)*realisedV;
                %EdU(ixt, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;
            end %ixY

        end %ixA
    end %ixAIME
    fprintf('Passed period %d of %d.\n',ixt, T)
end %ixt

%Check the stochastic discount factors
%checkSDF(policyA1, EV, dU);

end %function