% This file fits a determinatally-thinned point process to a
% (dependently-)thinned-point process based on the method outlined in the
% paper by Blaszczyszyn and Keeler[1], which is essentially the method
% developed by Kulesza and Taskar[2] in Section 4.1.1.
%
% This is the second file (of three files) to run to reproduce the results
% presented in the paper by Blaszczyszyn and Keeler[1].
%
% The data used for fitting (or training) is stored in the file Subset.mat,
% which is generated with the MATLAB file SubsetGenerate.m.
%
% The fitting paramters are stored locally in the file SubsetFitParam.mat
%
% REQUIREMENTS:
% Uses Optimization Toolbox (uses  MATLAB function fminsearch).
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
% [2] Kulesza and Taskar, "Determinantal point processes for machine
% learning",Now Publisers, 2012
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

close all
clearvars;
clc;

load('Subset.mat');

T=100;%Number of training/learning samples
options=optimset('Display','iter'); %options for fminsearch

%%%START -- Model fitting parameters -- START
%The parameters can be changed, but these were observed to work well.
if any(choiceModel==1:2)
    %Fitting parameters for Matern hard-core (I/II) point process
    N=1; %number of neighbours for distances -- use N=1 or N=2
    M=0;%Must be N=M-1 or M=0.
    booleOptSigma=true; %Set to true to also optimize sigma value
    choiceKernel=1; %1 for Gaussian kernel, 2 for Cauchy kernel
    sigma=1;%sigma for Gaussian or Cauchy kernel
    %Sigma value is ignored if booleOptSigma=true
else
    % Fitting parameters for triangle point process
    N=2; %number of neighbours for distances -- use N=1 or N=2
    M=1;%Must be N=M-1 or M=0.
    booleOptSigma=false; %Set to true to also optimize sigma value
    choiceKernel=1; %1 for Gaussian kernel, 2 for Cauchy kernel
    sigma=0;%sigma for Gaussian or Cauchy kernel
    %Sigma value is ignored if booleOptSigma=true
end
%%%END -- Model fitting parameters -- END

%Probability of a Poisson realization with two few points
probLessM=sum(poisspdf(0:N,lambda));
if any([ppStructPoisson(:).n]<=N)
    error('Underlying Poisson realization needs at least N points');
end
%total number of possible training/learning samples
numbTrain=size(ppStructPoisson,1);
if T>numbTrain
    error('Not enough training samples ie T>numbSim.');
end

%Deterministic (ie gradient) optimization method
if booleOptSigma
    thetaGuess=ones(N+M+2,1); %Take initial guess for sigma values
else
    thetaGuess=ones(N+M+1,1); %Take initial guess for theta values
end
%function to maximize. See below for funLikelihood_Data function
funMax_theta=@(theta)funLikelihood_data(T,ppStructPoisson,...
    indexCellSub,choiceKernel,lambda,sigma,theta,booleOptSigma,N,M);
%define function to be minimized
funMin=@(theta)(-1*funMax_theta(theta));
%Minimize function -- may take a while.
thetaMax=fminunc(funMin,thetaGuess,options); %minimize function

if booleOptSigma
    sigma=thetaMax(end); %retrive sigma values from theta
    thetaMax=thetaMax(1:end-1);
end
sigma
thetaMax

choiceModelFitted=choiceModel; %record which model was used for fitting

save('SubsetFitParam.mat','thetaMax','T','sigma','N','M', ...
    'choiceModelFitted','booleOptSigma','choiceKernel');

%Function definitions for log-likelihood.
function logLikelihood=funLikelihood_data(T,ppStructPoisson, ...
    indexCellSub,choiceKernel,lambda,sigma,theta,booleOptSigma,N,M)
if booleOptSigma
    %sets sigma to one of the theta parameters to optimize
    sigma=theta(end);
    theta=theta(1:end-1);
end
%initialize  vector
logLikelihoodVector=zeros(T,1);

%Loop through all training/learning samples
for tt=1:T
    xx=ppStructPoisson(tt).x;yy=ppStructPoisson(tt).y;
    indexSub=indexCellSub{tt}; %index for sub point process
    
    %Create L matrix (ie for Phi) based on nearest neighbours
    L=funNeighbourL(xx,yy,lambda,choiceKernel,sigma,theta,N,M);
    
    %Create sub L matrix (ie for Psi)
    subL=L(indexSub,indexSub);
    logLikelihoodVector(tt)=(log(det(subL))-log(det(L+eye(size(L)))));
end
logLikelihood=sum(logLikelihoodVector);
end
