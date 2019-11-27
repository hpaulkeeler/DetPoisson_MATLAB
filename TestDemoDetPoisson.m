% Randomly simulates a determinantally-thinned Poisson point process. It
% then tests the empirical results against analytic ones.
%
% A determinantally-thinned Poisson point process is essentially a discrete
% determinantal point process whose underlying state space is a single
% realization of a Poisson
% point process defined on some bounded continuous space.
%
% For more details, see the paper by Blaszczyszyn and Keeler[1].
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018.
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

close all;
clearvars; clc;

%set random seed for reproducibility
rng(1);

numbSim=10^5; %number of simulations

%%%START -- Parameters -- %%%START
%Poisson point process parameters
lambda=10; %intensity (ie mean density) of the Poisson process

%choose kernel
choiceKernel=1; %1 for Gaussian (ie squared exponetial );2 for Cauchy
sigma=1;% parameter for Gaussian and Cauchy kernel
alpha=1;% parameter for Cauchy kernel

%Simulation window parameters
xMin=0;xMax=1;yMin=0;yMax=1;
xDelta=xMax-xMin;yDelta=yMax-yMin; %rectangle dimensions
areaTotal=xDelta*yDelta; %area of rectangle
%%%END -- Parameters -- %%%END

%%% START -- Simulate a Poisson point process on a rectangle %%%START
%Simulate Poisson point process
%numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
numbPoints=5;
xx=xDelta*(rand(numbPoints,1))+xMin;%x coordinates of Poisson points
yy=xDelta*(rand(numbPoints,1))+yMin;%y coordinates of Poisson points
%%% END -- Simulate a Poisson point process on a rectangle --END %%%

%%% START -- CREATE L matrix -- START %%%
%all squared distances of x/y difference pairs
xxDiff=bsxfun(@minus,xx,xx');
yyDiff=bsxfun(@minus,yy,yy');
rrDiffSquared=(xxDiff.^2+yyDiff.^2);
if choiceKernel==1
    %%Gaussian/squared exponential kernel
    L=lambda*exp(-(rrDiffSquared)/sigma^2);
elseif choiceKernel==2
    %%Cauchy kernel
    L=lambda./(1+rrDiffSquared/sigma^2).^(alpha+1/2);
else
    error('choiceKernel has to be equal to 1 or 2.');
end
%%% END-- CREATE L matrix -- %%% END

%%% START Testing DPP simulation START %%%
%Retrieve eigenvalues and eigenvectors
[eigenVectL,eigenValL]=eig(L); %eigen decomposition

%run simulations with tests
probX_i_Emp=zeros(numbPoints,1); %initialize variables
indexTest=2:3; %choose a subset of [1 numbPoints]
probTestEmp=0; %initialize variables
%loop through for each simulation
for ss=1:numbSim
    %run determinantal simuation
    indexDPP=funSimSimpleLDPP(eigenVectL,eigenValL); %returns index
    probX_i_Emp(indexDPP)=probX_i_Emp(indexDPP)+1;

    countTemp=0; %initialize count
    for ii=1:length(indexTest)
        %check that each point of test subset appears
        countTemp=countTemp+any(indexDPP==indexTest(ii));
    end
    probTestEmp=probTestEmp+(countTemp==length(indexTest));
end

%empirically estimate the probabilities of each point appearing
probX_i_Emp=probX_i_Emp/numbSim

%calculate exactly the probabilities of each point appearing
K=funLtoK(L);
probX_i_Exact=diag(K)

%empirically estimate the probabilities of test subset appearing
probTestEmp=probTestEmp/numbSim

%calculate exactly the probabilities of test subset appearing
probTestExact=det(K(indexTest,indexTest))

%%%END Testing DPP simulation END%%%
