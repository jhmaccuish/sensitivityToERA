function [ y, c, a, v, l ] = simWithUncer(policyA1,policyL,EV,startingA)

% This function takes the policy functions and value functions, along with
% starting assets and returns simulated paths of income, consumption,
% assets and value

%% ------------------------------------------------------------------------ 
% Declare global we need this file have access to
global mu sigma rho T r Tretire
global Agrid  Ygrid numSims 
global normBnd benefit delta spouseInc numAIME pension

%% ------------------------------------------------------------------------
% Initialise arrays that will hold the paths of income consumption, value
% and assets

% Arguments for output
y = NaN(T, numSims);            % income
c = NaN(T, numSims);            % consumption
l = NaN(T, numSims);            % labour supply
v = NaN(T, numSims);            % value
a = NaN(T + 1,numSims);         % this is the path at the start of each period, so we include the 'start' of death

% Other arrays that will be used below
e = NaN(T, numSims);            % the innovations to log income
logy1 = NaN(1, numSims);        % draws for initial income
ly = NaN(T, numSims);           % log income
ypathIndex = NaN(T, numSims);   % holds the index (location) in the vector 

%% ------------------------------------------------------------------------
% Obtain time series of incomes for our simulated individuals
%-------------------------------------------------------------------------%#
    % Draw random draws for starting income and for innovations
    seed1 = 1223424; % For the innovations
    seed2 = 234636;  % For initial income
    sig_inc = sigma/ ((1-rho^2)^0.5);
    sig_initial = mean([0.073 0.053 0.110 0.112])^0.5;
    [ e ] = getNormalDraws( 0, sigma,  T, numSims, seed1);  % normally distributed random draws for the innovation
    [ logy1 ] =  getNormalDraws( mu, sig_initial,  1, numSims, seed2); % a random draw for the initial income       
  
    % Get all the incomes, recursively
    for s = 1:1: numSims                           % loop through individuals
         ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc );
         y(1, s) = exp(ly(1, s)+polyval(delta,1))+ spouseInc;      
          %y(1, s) = (s==1)*min(Ygrid(1,:))+(s>1)*max(Ygrid(1,:));
        for t = 1:1:T                              % loop through time periods for a particular individual               
            if (t ~= T)  % Get next year's income
                %y(t+1, s) = (s==1)*min(Ygrid(t,:))+(s>1)*max(Ygrid(t,:));
                ly(t+1, s) = (1 -rho) * mu + rho * ly(t, s) + e(t + 1, s);
                ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc );
                y(t+1, s) = exp( ly(t+1, s) + polyval(delta,t+1) );                
            end % if (t ~= T)
            y(t, s) = y(t, s)+(t >= Tretire)*pension+(t >= 45)*pension+ spouseInc;
       end % t     
    end % s

%% ------------------------------------------------------------------------
% Obtain consumption, asset and value profiles
%-------------------------------------------------------------------------%#  
     disp('Iteration:          \n');
     for s = 1:1: numSims
        fprintf('\b\b\b\b\b%5.0f',s);
        a(1, s) = startingA;                   
         for t = 1:1:T                              % loop through time periods for a particular individual               
            clear tA1 tV;                      %necessary as the dimensions of these change as we wor through this file
            tA1(:,:) = policyA1(t, :, :);  % the relevant part of the policy function
            tV(:,:) = EV(t, :, :);         % the relevant part of the value function
            %if (t < Tretire)                       % first for the before retirement periods
%             else                          % next for the post retirement periods     
%                 a(t+1, s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tA1, a(t, s), y(t, s));                    
%                 v(t  , s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tV , a(t, s), y(t, s));
% %                 a(t+1, s) = interp1(Agrid(t,:)', tA1, a(t, s), 'linear', 'extrap');
% %                 v(t,   s) = interp1(Agrid(t,:)', tV , a(t, s), 'linear', 'extrap');
%             end 
            [~,idxY]=min(abs(y(t, s)-Ygrid(t,:)));
            [~,idxA]=min(abs(a(t, s)-Agrid(t,:)));
            l(t,s)=policyL(t,idxA,idxY);
            if ~l(t,s)
                y(t, s)=benefit(t);
            end
            a(t+1, s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tA1, a(t, s), (y(t, s)));                    
            v(t  , s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tV , a(t, s), (y(t, s)));
            % Check whether next period's asset is below the lowest
            % permissable
            if ( a(t+1, s) < Agrid(t+1, 1))
               [ a(t+1, s) ] = checkSimExtrap( Agrid(t+1, 1),y(t, s), t ); 
            end
            
            % Get consumption from today's assets, today's income and
            % tomorrow's optimal assets
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+r));
        end   %t      
    end % s

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------% 
 
 
end