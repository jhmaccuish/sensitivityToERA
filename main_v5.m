%% ------------------------------------------------------------------------ 
% Dynamic Economics in Practice. 
% Monica Costa Dias and Cormac O'Dea
% Institute for Fiscal Studies
% 17th - 18th September 2014


%% ------------------------------------------------------------------------ 
% DESCRIPTION
% This program solves and simulates a finite period consumption and saving 
% problem. There is income that can be uncertain. The income program can be
% hardcoded by the user or can be set to follow a log normal autoregessive
% process

%% ------------------------------------------------------------------------ 
% PREAMBLE
% Ensure that all storage spaces variables, matrices, memory, globals are 
% clean from information stored in past runs

tic;        % start the clock
%clear;  % clear memory
close all;  % close any graphs


%% ------------------------------------------------------------------------ 
% DECLARE VARIABLES AND MATRICES THAT WILL BE 'GLOBAL'
% explicitly set variables and matrices to be shared throughout the routine
% as globals

global beta gamma r T sigma mu rho Tretire          % structural model parameters (1)
global isUncertainty                                % structural model parameters (2)
global delta
global numPointsA Agrid                             % assets grid and dimension
global numPointsY Ygrid incTransitionMrx            % income grid, dimension, and transition matrices (1)
global interpMethod linearise                       % numerical methods to be used
global tol minCons numSims normBnd                  % numerical constants
global plotNumber                                   % for graphing
%Added
global numPointsL                                   % number of labour market states
global benefit                                      % benefit
global pension                                      % pension
global nu                                           % labour-leisure relative weight
global hrsWrk                                       % hours worked
global probU                                        % prob unemployment
global probE
global spouseInc
global numAIME AIMEgrid

%% ------------------------------------------------------------------------ 
% NUMERICAL METHODS
% select solution, interpolation and integration methods

interpMethod = 'pchip';      % interpolation method - choose from 'linear', 'nearest', 'spline', 'pchip'
solveUsingValueFunction = 1; % solution method: set to 1 to solve using value function, else =0
solveUsingEulerEquation = 0; % solution method: set to 1 to solve using Euler equation, else =0
linearise = 1;               % whether to linearise the slope of EV when using EE - set linearise=1 to do this, else = 0


%% ------------------------------------------------------------------------ 
% NUMERICAL CONSTANTS
% set constants needed in numerical solution and simulation

% precision parameters
%--------------------------------------%
tol = 1e-10;                 % max allowed error
minCons = 1e-5;              % min allowed consumption

% where to truncate the normal distributions
%--------------------------------------%
normBnd = 4;                 %Ignore draws less than -NormalTunc*sigma and greater than normalTrunc*sigma 

% information for simulations
%--------------------------------------%
numSims = 10000;                % How many individuals to simulate


%% ------------------------------------------------------------------------ 
% THE ECONOMIC ENVIRONMENT
% Set values of structural economic parameters
T = 75;                      % Number of time period
r = 0.01;                    % Interest rate
beta = 0.970;                % Discount factor
gamma = 3;%1.5;              % Coefficient of relative risk aversion
mu = 0;                      % mean of initial log income
sigma = mean([0.007 0.005 0.010 0.012])^0.5;%0.073;                % variance of innovations to log income
rho = 0.96;%0.75;                  % persistency of log income
Tretire = 41;                % age after which there is no income earned
borrowingAllowed = 0;        % Is borrowing allowed
isUncertainty = 1;           % Is income uncertain?
startA = 0;                  % How much asset do people start life with
%Added parameters
nu = mean([0.422,0.417,0.5167,0.499]); %lesuire
hrsWrk = 1840/(5824);
delta(3) = mean([9.083298	8.5256042	7.906901	7.392151]);%0;%9.083298;%
delta(2)=mean([0.029333	0.061582208	0.094386	0.121426]);%0;%0.029333;%
delta(1)=mean([-0.00023	-0.0005901	-0.00091	-0.00117]);%0;%-0.00023;%
probU = mean([0.227 0.099 0.201 0.063;
              0.049 0.026 0.024 0.016;
              0.027 0.010 0.054 0.01;
              0.014 0.013 0.032 0.006;
              0.004 0.008 0.029 0.003;
              0.007 0.012 0.038 0.002;
              0.013 0.003 0.020 0.004;
              0.017 0.003 0.005 0.001;
              0.007 0.002 0.004 0.003;
              0.010 0.000 0.010 0.000;],2);
probE = mean([0.130 0.174 0.118 0.238;
              0.043 0.005 0.029 0.032;
              0.022 0.013 0.037 0.020;
              0.005 0.017 0.011 0.008;
              0.002 0.013 0.008 0.008;
              0.007 0.017 0.013 0.012;
              0.010 0.013 0.003 0.008;
              0.002 0.003 0.000 0.004;
              0.000 0.003 0.005 0.000;
              0.007 0.000 0.000 0.004],2);
spouseInc = mean([5121, 5282, 6840, 7684]);          
%delta(3)=delta(3)+spouseInc;
probE = kron(probE, [1;1]);
probU = kron(probU, [1;1]);

%% ------------------------------------------------------------------------ 
% GRIDS
% choose dimensions, set matrices and select methods to construct grids

%The grid for assets
%--------------------------------------%
numPointsA = 20;             % number of points in the discretised asset grid
gridMethod = '3logsteps';    % method to construct grid. One of equalsteps, logsteps, 3logsteps, 5logsteps or 10logsteps

%The grid for income shocks
%--------------------------------------%
numPointsY = 20;           %  points in grid for income (should be 2 if hard-coded)
numAIME = 10;

%The grid labour choices
%--------------------------------------%
numPointsL = 2;           %  points in grid for income (should be 2 if hard-coded)

%% Check inputs
checkInputs;

%% Get income grid
%pension =   2.5491e+03;%3.1425e+03;%0.2649;%0.1968;
[Ygrid, incTransitionMrx, minInc, maxInc, AIMEgrid] = getIncomeGrid; 
%benefit = benefit;
%pension = 0.17*mean(mean(Ygrid));

%% ------------------------------------------------------------------------ 
% GET ASSET GRID
% populate grid for assets using 'gridMethod'
% populate matrix borrowingConstraints


[ borrowCon, maxAss ] = getMinAndMaxAss(borrowingAllowed, benefit, maxInc, startA);

Agrid = NaN(T+1, numPointsA);
for ixt = 1:1:T+1
    Agrid(ixt, :) = getGrid(borrowCon(ixt), maxAss(ixt), numPointsA, gridMethod);
end


%% ------------------------------------------------------------------------ 
% SOLVE CONSUMER'S PROBLEM
% Get policy function and value function 

if solveUsingValueFunction == 1
    [ policyA1, policyC, policyL, val, exVal, EdU ] = solveValueFunction;
elseif solveUsingEulerEquation == 1
    %[ policyA1, policyC, val, exVal, ] = solveEulerEquation;
end
% toc;
% return;
%% ------------------------------------------------------------------------ 
% SIMULATE CONSUMER'S PATHS
% start from initial level of assets and simulate optimal consumption and
% savings profiles over lifecycle

fprintf('Start simulation\n')
if isUncertainty == 0
    %[ ypath, cpath, apath, vpath ] = simNoUncer(policyA1, exVal, startA);
else
    [ ypath, cpath, apath, vpath, lpath ] = simWithUncer(policyA1,policyL,exVal, startA);
end

% %  figure();
% %  plot(21:(20+T),mean(lpath,2));
% %  title('Mean Participation');
%
% %[~,order]=sort(mean(apath));
% 
 %[~,order]=sort(mean(apath,1));
 [~,order]=sort(apath(40,:));
 oLpath=lpath(:,order);
 oApath=apath(:,order);
 oCpath=cpath(:,order);
 
 plot(21:(20+T),mean(oLpath(:,end-500:end),2));
 title('Mean Participation Richest 10% by life-time wealth');
toc;
% same=zeros(1,T);
% for i=1:T;
%   same(i)=min(min(policyLCopy65(i,:,:)==policyL(i,:,:))) ; 
% end;
% sameA=zeros(1,T);
% for i=1:T;
%   sameA(i)=min(min(policyA1Copy65(i,:,:)==policyA1(i,:,:))) ; 
% end;

% toc;
% return;
% %% ------------------------------------------------------------------------ 
% % PLOTS
% % Plots some features of the simulation and simulation
% 
% %Plot paths of consumption, income and assets 
% plotNumber = 0;
% plotCpath(cpath)
% plotApath(apath, borrowCon)
% plotYAndCpaths( ypath, cpath );
% plotYCAndApaths( ypath, cpath, apath );
% 
% % Now plot value and policy functions
% whichyear = 20;
% plotNode1 = 3;
% plotNodeLast = numPointsA; 
% plots; 
% % 
% % toc;     % Stop the clock
% % % ------------------------------------------------------------------------ 
% ------------------------------------------------------------------------   