% Randomly simulates a determinantally-thinned Poisson point process.
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

%%%START -- Parameters -- %%%START
%Poisson point process parameters
lambda=50; %intensity (ie mean density) of the Poisson process

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
numbPoints=poissrnd(areaTotal*lambda);%Poisson number of points
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

%START - Sampling/simulating DPP - START
%Retrieve eigenvalues and eigenvectors
[eigenVectorsL,eigenValuesL]=eig(L); %eigen decomposition
eigenVectorsK = (eigenVectorsL);
eigenValuesL=abs(diag(eigenValuesL)); %eigenvalues of L
eigenValuesK = eigenValuesL./(1+eigenValuesL); %eigenvalues of K

%Bernoulli trials (ie coin flips) to determine number of points
indexEigen = find(rand(length(eigenValuesK),1) <= eigenValuesK);

%number of points in the DPP realization
numbPointsDPP=length(indexEigen);
%retrieve eigenvectors corresponding to successful Bernoulli trials
spaceV = eigenVectorsK(:,indexEigen); %subspace V
indexDPP = zeros(numbPointsDPP,1); %index for final DPP configuration

%Loop through for all points
for ii=numbPointsDPP:-1:1
    %Compute probabilities for each point i
    Prob_i=sum(spaceV.^2,2); %sum across rows
    Prob_i=Prob_i/sum(Prob_i); %normalize

    %Choose a point (from 1 to numbPoints) using (prob mass function) Prob_i
    indexCurrent=find(rand(1)<=cumsum(Prob_i),1);
    indexDPP(ii)=indexCurrent;

    %Choose a vector to remove
    jj=find(spaceV(indexCurrent,:),1);
    columnVj=spaceV(:,jj); %j-th column of V
    spaceV=spaceV(:,[1:jj-1 jj+1:end]);

    %Update matrix V by removing Vj component from the space
    spaceV=spaceV-columnVj.*spaceV(indexCurrent,:)/columnVj(indexCurrent);

    %Orthonormalize using Householder method
    [spaceV,~]=qr(spaceV,0);
end
indexDPP=sort(indexDPP); %sort index in ascending order
%END - Sampling/simulating DPP - END

%%% START -- Plotting -- START %%%
figure;hold on;
markerSizeNumb=80; %marker size of markers colors
vectorColor=rand(1,3).^(1); %random vector for colors of
%Plot Poisson point process
plot(xx,yy,'ko','MarkerSize',markerSizeNumb/6);
%Plot determinantally-thinned Poisson point process
plot(xx(indexDPP),yy(indexDPP),'.','MarkerSize',markerSizeNumb/3,'color',vectorColor);
grid;
axis square;set(gca,'YTick',[]); set(gca,'XTick',[]);
legend('Poisson process', 'Determinantal Poisson');
%%% END -- Plotting -- END %%%
