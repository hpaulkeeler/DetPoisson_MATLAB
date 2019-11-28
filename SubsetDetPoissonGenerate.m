% This simulates determinatally-thinned point processes that have been 
% fitted to thinned-point process based on the method outlined in the paper 
% by Blaszczyszyn and Keeler[1], which is essentially the method developed 
% by Kulesza and Taskar[2]. It then gathers empirical estimates of nearest
% neighbour distributions and contact distributions; for details see,
% the books by Chiu, Stoyan, Kendall, and Mecke[3] or Baddeley, Rubak and
% Turner[4]. It also calculates these distributions using determinant
% equations by simulating the underlying Poisson point process, and not the
% determinatally-thinned point process; see Blaszczyszyn and Keeler[1],
% Section III. Part of this is based on the (reduced) Palm distribution of
% determinantal point processes derived by Shirai and Takahashi[5]
%
% This is the third file (of three files) to run to reproduce the results
% presented in the paper by Blaszczyszyn and Keeler[1].
%
% The data used for fitting (or training) is stored in the file Subset.mat,
% which is generated with the MATLAB file SubsetGenerate.m.
%
% The fitting paramters are stored locally in the file SubsetFitParam.mat
%
% REQUIREMENTS:
% Uses Statistics (and Machine learning) Toolbox for random variable.
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
% [2] Kulesza and Taskar, "Determinantal point processes for machine
% learning",Now Publisers, 2012
% [3] Chiu, Stoyan, Kendall, and Mecke, "Stochastic geometry and its
% applications", Wiley.
% [4] Baddeley, Rubak and Turner, "Spatial point patterns: Methodology and
% applications with R, 2016.
% [5] Shirai and Takahashi, "Random point fields associated with certain
% Fredholm determinants I -- fermion, poisson and boson point", 2003.
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

close all
clearvars;
clc;

%seedRand=1; rng(seedRand); %Seed for random number generation

boolePlot=true; %plot configurations
numbSim=10^3;
booleEquation=true; %if booleEquation=true, calculate contact and nearest
%neighbour distributions using equation, and not just empirical estimate.

%Load data from .mat files created by other .m scripts
load('Subset.mat'); %Created by SubsetGenerate.m
load('SubsetFitParam.mat'); %Created by SubsetDetPoissonFit.m

%Helper function: distance (squared) between points
funDistSquared=@(x,y,u,v)((x-u).^2+(y-v).^2);

n=(size(ppStructPoisson,1));%total number of training/learning samples
if (numbSim+T>n)
    error('Need to create more realizations with SubsetGenerate.m');
else
    %Look at unused realizations (ie the ones not used for fitting)
    %select a random subset of unused realizations
    ttValuesPerm=(randperm(numbSim)+T);
    ttValues=ttValuesPerm(1:numbSim);
end

if booleDisk %disk simulation/sample/window
    xxCenter=xx0;  yyCenter=yy0; %center of sample window
    rMax=windowSample;
else %rectangle
    %center of sample window
    xxCenter=(windowSample(2)-windowSample(1))/2;
    yyCenter=(windowSample(4)-windowSample(3))/2;
    %find largest radius for bounding disk
    rMax=min(windowSample(2)-windowSample(1),...
        windowSample(4)-windowSample(3));
end

%initialize  variables for collecting statistics
meanNumbDPPCond=zeros(numbSim,1);
numbSub=zeros(numbSim,1);
numbDPP=zeros(numbSim,1);
distContactSub=Inf(numbSim,1);distContactDPP=Inf(numbSim,1);
distNearestMidSub=Inf(numbSim,1);distNearestMidDPP=Inf(numbSim,1);

if booleEquation
    %ininiate variables if Palm probability expressions are used
    rNearestDPP_Eq=linspace(0,rMax,30);
    GDPP_Eq=zeros(size(rNearestDPP_Eq));
    rContactDPP_Eq=linspace(0,rMax,30);
    HDPP_Eq=zeros(size(rContactDPP_Eq));
    pi_x=0;
end
%create cell array for all nearest neighbour distances
distNearestAllSubCell=cell(numbSim,1);
distNearestAllDPPCell=cell(numbSim,1);
for ss=1:numbSim
    tt=ttValues(ss);
    xxPoisson=ppStructPoisson(tt).x;yyPoisson=ppStructPoisson(tt).y;
    numbRandPoisson=ppStructPoisson(tt).n;
    indexSub=indexCellSub{tt}; %index for sub point process
    xxSub=xxPoisson(indexSub);yySub=yyPoisson(indexSub);
    numbSub(ss)=length(indexSub);
    %generate L matrix based on Poisson point realization
    L=funNeighbourL(xxPoisson,yyPoisson,lambda,...
        choiceKernel,sigma,thetaMax,N,M);
    %eigen decomposition of L
    [eigenVectorsL,eigenValuesL]=eig(L);
    eigenVL=diag(eigenValuesL);
    eigenValuesK=eigenVL./(1+eigenVL);%find K eigenvalues
    %average number of DPP points
    meanNumbDPPCond(ss)=mean(real(eigenValuesK))*numbRandPoisson;
    %Simulate next DPP generation
    indexDPP=funSimSimpleLDPP(eigenVectorsL,eigenValuesL);
    %assign points
    xxDPP=xxPoisson(indexDPP);yyDPP=yyPoisson(indexDPP);
    numbDPP(ss)=length(xxDPP);%number of points
    
    %%% START -- Empirical distributions of contact and nearest-neighbour
    %Nearest neighbour chooses a point closest to the centre, which
    %introduces a small bias.
    %subset (ie Matern) point process
    if (~isempty(indexSub)>0)
        [distContactSub(ss),indexCentreSub]=min(funDistSquared(xxSub,yySub,xxCenter,yyCenter));
        %nearest neighbour stats (only collect if there is more than one point)
        if (length(indexSub)>1)
            %subset (ie Matern) point process
            indexTemp=setdiff(1:numbSub(ss),indexCentreSub);
            distNearestMidSub(ss)=...
                min(funDistSquared(xxSub(indexTemp),yySub(indexTemp)...
                ,xxSub(indexCentreSub),yySub(indexCentreSub)));
            %all nearest neighbour distances
            [~,distNearestSubTemp]=knnsearch([xxSub,...
                yySub],[xxSub,yySub]...
                ,'K',2);
            distNearestAllSubCell(ss)={distNearestSubTemp(:,2)};
        end
    end
    %DPP process
    if (~isempty(indexDPP)>0)
        [distContactDPP(ss),indexCentreDPP]=min(...
            funDistSquared(xxDPP,yyDPP,xxCenter,yyCenter));
        %nearest neighbour stats (only collect if there is more than one point)
        if (length(indexDPP)>1)
            indexTemp=setdiff(1:numbDPP(ss),indexCentreDPP);
            distNearestMidDPP(ss)=...
                min(funDistSquared(xxDPP(indexTemp),yyDPP(indexTemp)...
                ,xxDPP(indexCentreDPP),yyDPP(indexCentreDPP)));
            %all nearest neighbour distances
            [~,distNearestDPPTemp]=knnsearch([xxDPP,...
                yyDPP],[xxDPP,yyDPP]...
                ,'K',2);
            distNearestAllDPPCell(ss)={distNearestDPPTemp(:,2)};
        end
    end
    if booleEquation
        %Calculate Palm distriubtion for a point U existing at (xxU,yyU)
        xxU=0;yyU=0;
        %add point to Poisson realization
        xxPoisson0=[xxU;xxPoisson]; yyPoisson0=[yyU;yyPoisson];
        indexPalm=1;%index for point U
        %calculate Palm version of L matrix
        L0=funNeighbourL(xxPoisson0,yyPoisson0,lambda,...
            choiceKernel,sigma,thetaMax,N,M);%find Palm L matrix
        %convert L into K matrix
        K=L/(eye(size(L))+L);
        K0=L0/(eye(size(L0))+L0);
        pi_x=pi_x+K0(indexPalm)/numbSim; %estimate thinnning density
        %Palm (ie Shira and Taskahsi) construction
        KPalm=K0-K0(:,indexPalm).*K0(indexPalm,:)./K0(indexPalm,indexPalm);
        KPalmReduced=KPalm(2:end,2:end); %reduced Palm version
        for jj=1:length(rNearestDPP_Eq)
            r=rNearestDPP_Eq(jj);
            booleBall=funDistSquared(xxPoisson,yyPoisson,xxU,yyU)<r^2;
            GDPP_Eq(jj)=GDPP_Eq(jj)+(1-det(eye(sum(booleBall)) ...
                -KPalmReduced(booleBall,booleBall)))*K0(indexPalm);
            r=rContactDPP_Eq(jj);
            booleBall=funDistSquared(xxPoisson,yyPoisson, ...
                xxCenter,yyCenter)<r^2;
            HDPP_Eq(jj)=HDPP_Eq(jj)+det(eye(sum(booleBall))-...
                K(booleBall,booleBall));
        end
        
    end
    %%%END -- Empirical distributions of contact and nearest-neighbour
    
end
distNearestAllSub=cell2mat(distNearestAllSubCell); %convert to vector
distNearestAllDPP=cell2mat(distNearestAllDPPCell); %convert to vector

%square root the squared distances
distContactSub=sqrt(distContactSub);
distContactDPP=sqrt(distContactDPP);
distNearestMidSub=sqrt(distNearestMidSub);
distNearestMidDPP=sqrt(distNearestMidDPP);

%calculate empirical distributions
[HSub,rContactSub]=ecdf(distContactSub);
[HDPP,rContactDPP]=ecdf(distContactDPP);
[GMidSub,rNearestMidSub]=ecdf(distNearestMidSub);
[GMidDPP,rNearestMidDPP]=ecdf(distNearestMidDPP);
[GAllSub,rNearestAllSub]=ecdf(distNearestAllSub);
[GAllDPP,rNearestAllDPP]=ecdf(distNearestAllDPP);

%parameters
lambda
sigma
thetaMax

%Empirical density
lambdaNumDPP=mean(meanNumbDPPCond)/areaSample
lambdaEmpDPP=mean(numbDPP)/areaSample
lambdaEmpSub=mean(numbSub)/areaSample

%average nearest neighbour
meanDistNearestMidSub=mean(distNearestMidSub)
meanDistNearestMidDPP=mean(distNearestMidDPP)
meanDistNearestAllSub=mean(distNearestAllSub)
meanDistNearestAllDPP=mean(distNearestAllDPP)

%%%START -- PLOTTING -- START%%%
vectorColor1=[1,0,0];%rand(1,3);%random colour vector
vectorColor2=[0,0,1];%rand(1,3).^2;%random colour vector

labelModelPlot=strcat(labelModel,' (empirical)');
figure;
%plot empirical distribution of nearest neighbour distance

set(gcf,'DefaultLineLineWidth',1.2);
plot(rNearestAllSub,GAllSub,'-','Color',vectorColor1);
hold on;
plot(rNearestAllDPP,GAllDPP,'--','Color',vectorColor2);
grid;
if booleEquation
    %calculate nearest neighbour distribution based on Palm probability
    GDPP_Eq=GDPP_Eq(:)/numbSim/pi_x;
    meanDistNearestDPP_eq=trapz(rNearestDPP_Eq(:).*[0;diff(GDPP_Eq)])
    plot(rNearestDPP_Eq,GDPP_Eq,'+-','Color',vectorColor2);
    legend(labelModelPlot,...
        'Determinatally-thinned (empirical)',...
        'Determinatally-thinned (equation)',...
        'Location','SouthEast','AutoUpdate','off');
else
    legend(labelModelPlot,...
        'Determinatally-thinned (empirical)',...
        'Location','SouthEast');
end
if (choiceModel<=2)
    plot(rSub,0,'r'); %indicate Matern radius
end
axis([0,rMax,0,1]);
xlabel('r','fontsize',16);  ylabel('G(r)','fontsize',16);
print('-clipboard','-dbitmap'); %copy to clipboard

if (boolePlot==1)
    %%% START Plotting %%%
    %Plot intial results
    markerSizeNumb=50; %marker size of markers colors
    markerSizeRatio=1;
    figure;
    set(gcf,'DefaultLineLineWidth',1.2);
    %draw underlying Poisson point process
    plot(xxPoisson,yyPoisson,'ko','MarkerSize',markerSizeNumb/4.0);hold on;
    %example of a subset (eg Matern hard-core) process
    plot(xxSub,yySub,'x','MarkerSize',markerSizeNumb/6.0,'color',vectorColor1);
    %draw determinantal point process
    plot(xxDPP,yyDPP,'+','MarkerSize',markerSizeNumb/4,'color',vectorColor2);
    axis square;set(gca,'YTick',[]); set(gca,'XTick',[]);
    legend('Poisson process', labelModel, 'Determinantal Poisson');
    grid;
    
    %plot empirical distribution of contact distance
    figure;
    set(gcf,'DefaultLineLineWidth',1.2);
    plot(rContactSub,HSub,'-','Color',vectorColor1);
    hold on;
    plot(rContactDPP,HDPP,'--','Color',vectorColor2);
    grid;
    
    if booleEquation
        HDPP_Eq=1-HDPP_Eq(:)/numbSim;
        plot(rContactDPP_Eq,HDPP_Eq,'+-','Color',vectorColor2);
        
        legend(labelModelPlot,...
            'Determinatally-thinned (empirical)',...
            'Determinatally-thinned (equation)',...
            'Location','SouthEast');
    else
        legend(labelModelPlot,...
            'Determinatally-thinned (empirical)',...
            'Location','SouthEast');
    end
    xlabel('r','fontsize',16);  ylabel('H(r)','fontsize',16);
    axis([0,rMax,0,1]);
end
%%%END -- PLOTTING -- END%%%

%%% Disply some additional information
if choiceModelFitted==1
    disp('The fitted model is based on a Matern I point process.');
elseif choiceModelFitted==2
    disp('The fitted model is based on a Matern II point process.');
else
    disp('The fitted model is based on a Triangle point process.');
end
if choiceModel==1
    disp('The new test results are from a Matern I point process.');
elseif choiceModelFitted==2
    disp('The new test results are from a Matern II point process.');
else
    disp('The new test results are from a Triangle point process.');
end
if sigma~=0
    if choiceKernel==1
        disp('A Gaussian kernel was used for the similarity matrix S.');
    else
        disp('A Cauchy kernel was used for the similarity matrix S.');
    end
end
