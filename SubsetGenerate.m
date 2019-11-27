% This file can simulate three different types of dependently-thinned point
% processes based on an underling homogeneous Poisson point process. The 
% resulting realizations of these point processes are considered as subsets 
% of the underlying Poisson realizations. The results are stored locally in
% the file Subset.mat. It also plots a single realization.
%
% The two of the point processes are Matern I/II hard-core point
% processes. The other is a point process, which we call the triangle point
% process. All three are constructed by dependently thinning a homogeneous 
% Poisson point process.
% 
% It creates 10^4 realizations of Matern I/II point processes in a few
% seconds on a standard machine, whereas the triangle point process is 
% slower to generate.
%
% This is the first file (of three files) to run to reproduce the results
% presented in the paper by Blaszczyszyn and Keeler[1].
%
% The underlying Poisson realizations are stored as a structure (array)
% called ppStructPoisson, where the thinned point processes (represented as
% indices) are stored in a cell (array) called indexCellSub.
%
% For more details, see the paper by Blaszczyszyn and Keeler[1].
% For details on Matern hard-core point point processes, see, for example,
% books[2](pages 176--178), [3](pages 93--94) or [4](pages 37--38).
%
% REQUIREMENTS:
% Uses Statistics (and Machine learning) Toolbox for Poisson variable 
% generation (uses MATLAB function poissrnd).
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018.
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
% [2] Chiu, Stoyan, Kendall, and Mecke, "Stochastic geometry and its
% applications", Wiley.
% [3] Schneider and Weil, "Stochastic and integral geometry", Springer.
% [4] Blaszczyszyn, Haenggi, Keeler, and Mukherjee, "Stochastic geometry
% analysis of cellular networks", Cambridge University Press.
%
% Code available here:
% Keeler, 2018, https://github.com/hpaulkeeler/DetPoisson_MATLAB

close all; clearvars; clc;

%seedRand=1;rng(seedRand); %Seed for random number generation 
choiceModel=2;%1 for Matern I, 2 for Matern II; 3 for triangular process
booleDisk=true;%true (or 1) to sample points on a origin-centered disk
booleCheck=true; %change to true (or 1) to see empirical stations and plots
numbSim=10^4; %number of realizations/simulations created

%Point process parameters
lambda=10; %intensity (ie mean density) of underlying Poisson point process
rMatern=.8/sqrt(lambda); %inhibition radius of Matern I/II process
rTri=2/sqrt(lambda); %thinning parameter for "triangle" process
%thin points if triangle perimeter is smaller than 2*rTri

%sample window dimensions
xx0=0;yy0=0;
%disk centered at (xx0,yy0);
rDisk=1; %radius
%rectangle with bottom left corner at (xx0,yy0)
widthRect=1; %width
heightRect=1; %height
xMin=xx0;xMax=xMin+widthRect;yMin=yy0;yMax=yMin+heightRect;

booleExtended=true; %true to simulate Poisson points on extended window
constN=2;constM=constN-1; %additonal triangle process paramters

%Helper function: distance (squared) between points
funDistSquared=@(x,y,u,v)((x-u).^2+(y-v).^2);

if any(choiceModel==1:2)
    rExtended=rMatern;
elseif choiceModel==3
    rExtended=2*rTri;
else
    error('The variable choiceModel must be one of the following numbers: 1, 2 or 3');
end

if booleDisk %disk
    %Helper functions
    funArea=@(w)(pi*w.^2); %area function
    %check if points in window (returns boolean variables)
    funInWindow=@(x,y)(((x-xx0).^2+(y-yy0).^2)<rDisk^2);
    %polar to Cartesian
    funPol2Cart=@(theta,r)([r.*cos(theta)+xx0,r.*sin(theta)+yy0]);
    %Points uniformally colocated in a origin centred disk
    funUniPoints=@(n,w)funPol2Cart(2*pi*rand(n,1),w*sqrt(rand(n,1)));
    %create sample and big window
    windowSample=rDisk;
    windowExtended=rDisk+rExtended;
else %rectangle
    funArea=@(w)(abs(w(2)-w(1)).*abs(w(4)-w(3))); %area function
    %check if points in window (returns boolean variables)
    funInWindow=@(x,y)((x>=xMin)&(x<=xMax)&(y>=yMin)&(y<=yMax));
    %Points uniformally colocated in a origin corned rectangle
    funUniPoints=@(n,w)([(w(2)-w(1))*rand(n,1)+w(1), ...
        (w(4)-w(3))*rand(n,1)+w(3)]);
    %create sample and big window
    windowSample=[xMin,xMax,yMin,yMax];
    windowExtended=[xMin-rExtended,xMax+rExtended,...
        yMin-rExtended,yMax+rExtended];
end
if booleExtended
    windowBig=windowExtended; %use extended window for Poisson points
else
    windowBig=windowSample; %use sample window for Poisson points
    funInWindow=@(x,y)((ones(size(x))));
end

%calculate areas
areaSample=funArea(windowSample);
areaBig=funArea(windowBig);

%%%START -- Simulation section -- START%%%
%Marked Poisson pont process on bigWindow
%Generate *all* ensembles at once
pointsNumberBig=poissrnd(areaBig*lambda,numbSim,1);%Poisson number of points
pointsNumberTotal=sum(pointsNumberBig);
%uniform x/y coordinates of Poisson points
xxyyTemp=funUniPoints(pointsNumberTotal,windowBig);
xxAllBig=xxyyTemp(:,1); yyAllBig=xxyyTemp(:,2); %extract x/ypoints from matrix
markAge=rand(size(xxAllBig));%iid uniform marks for Matern II
%convert Poisson point processes into cells
xxCellBig=mat2cell(xxAllBig,pointsNumberBig);
yyCellBig=mat2cell(yyAllBig,pointsNumberBig);
numberCell=mat2cell(pointsNumberBig,ones(numbSim,1)); %cell for number of points
markCell=mat2cell(markAge,pointsNumberBig);
%create structure to represent Poisson point process
ppStructPoissonBig=struct('x',xxCellBig,'y',yyCellBig,...
    'n',numberCell,'window',windowBig,'marks',markCell);
ppStructPoisson=ppStructPoissonBig;

%Counts for estimating densities of generated point processes
pointCountsI=zeros(numbSim,1);pointCountsII=zeros(numbSim,1);
pointCountsTri=zeros(numbSim,1); pointCounts=zeros(numbSim,1);

%create cells
xxCellSub=mat2cell(zeros(numbSim,1),ones(numbSim,1)); %for x coordinates
yyCellSub=mat2cell(zeros(numbSim,1),ones(numbSim,1)); %for y coordinates
%index of x/y points for subset point process
indexCellSub=mat2cell(zeros(numbSim,1),ones(numbSim,1));
tic;
%Loop for every Poisson realization
for jj=1:numbSim
    xxBig=ppStructPoissonBig(jj).x;yyBig=ppStructPoissonBig(jj).y;
    markAge=ppStructPoissonBig(jj).marks;
    booleWindow=funInWindow(xxBig,yyBig);
    indexWindow=find(booleWindow);
    xx=xxBig(indexWindow); yy=yyBig(indexWindow);
    numbPoints=length(indexWindow);
    
    %Update underlying Poisson point process
    ppStructPoisson(jj).x=xx; ppStructPoisson(jj).y=yy;
    ppStructPoisson(jj).n=numbPoints;
    ppStructPoisson(jj).window=windowBig; %using big window
    ppStructPoisson(jj).marks=markAge(indexWindow);
    
    if choiceModel<=2
        %Matern I/II
        removeBooleI=false(size(xx));%Index for removing points -- Matern I
        keepBooleII=false(size(xx));%Index for keeping points -- Matern II
        for ii=1:numbPoints
            distanceTemp=funDistSquared(xx(ii),yy(ii),xxBig,yyBig); %distances to other points
            inDiskBoole=(distanceTemp<rMatern^2)&(distanceTemp>0);
            %Matern I
            removeBooleI(ii)=any(inDiskBoole);
            %Matern II
            %keep the younger points
            keepBooleII(ii)=all(markAge(indexWindow(ii))<markAge(inDiskBoole));
            %Note: if markAge(inDiskBoole) is empty, keepBooleII(ii)=1.
        end
        %Remove/keep points to generate Matern hard-core processes
        if choiceModel==1
            %Matern I
            keepBooleI=~(removeBooleI);
            xxMaternI=xx(keepBooleI);yyMaternI=yy(keepBooleI);
            indexCellSub{jj}=find(keepBooleI); %subset index for Matern I
            pointCountsI(jj)=numel(xxMaternI);
            xxCellSub{jj}=xxMaternI;
            yyCellSub{jj}=yyMaternI;
            
        else
            %Matern II
            xxMaternII=xx(keepBooleII);yyMaternII=yy(keepBooleII);
            indexCellSub{jj}=find(keepBooleII); %subset index for Matern II
            pointCountsII(jj)=numel(xxMaternII);
            xxCellSub{jj}=xxMaternII;
            yyCellSub{jj}=yyMaternII;
        end
        
    else        
        %START -- thinning points based on triangular perimiter -- START
        numbPointsBig=length(xxBig);
        booleMatrix=logical(ones(numbPointsBig)-eye(numbPointsBig)); %removes distances to self
        indexMatrix=repmat((1:numbPointsBig),numbPointsBig,1);
        %find every pair combination
        xxMatrix=repmat(xxBig,1,numbPointsBig-1);yyMatrix=repmat(yyBig,1,numbPointsBig-1);
        indexTemp=flipud(reshape(indexMatrix(booleMatrix),numbPointsBig,numbPointsBig-1));
        %calculate nearest distance for all pairs
        [distNearAll,indexNear]=sort(sqrt((xxMatrix-xxMatrix(indexTemp)).^2+(yyMatrix-yyMatrix(indexTemp)).^2),2);
        distNear=distNearAll(:,1:constN); %nearest distances
        %Run for distances between nearest neighbours
        %find indices of nearest neighbours
        indexBetween=indexTemp(sub2ind([numbPointsBig,numbPointsBig-1],repmat((1:numbPointsBig)',1,numbPointsBig-1),indexNear));
        distBetween=sqrt((xxMatrix(indexBetween(:,1:constM))-xxMatrix(indexBetween(:,2:constM+1))).^2 ...
            +(yyMatrix(indexBetween(:,1:constM))-yyMatrix(indexBetween(:,2:constM+1))).^2);
        %total of three distances ie the perimiter distance of the triangle
        distTriAll=sum(distNear,2)+distBetween;
        distTri=distTriAll(booleWindow);
        keepBooleTri=(distTri>rTri); %keep points with large triangle distance
        xxTri=xx(keepBooleTri);yyTri=yy(keepBooleTri);
        %END -- thinning points based on triangular perimeter -- END
        
        %index of points
        indexCellSub{jj}=find(keepBooleTri); %subset index for triangular
        pointCountsTri(jj)=numel(xxTri);
        xxCellSub{jj}=xxTri;
        yyCellSub{jj}=yyTri;
        
    end
    %%%END -- Simulation section -- END%%%
    
    %Collect statistics on number of points
    pointCounts(jj)=numel(indexWindow);
end
toc;

%%%START -- Create point process object -- START%%%
if (choiceModel==1)
    %area of hard-core disk
    diskMaternArea=pi*rMatern^2;
    %Matern I
    numberCell=mat2cell(pointCountsI,ones(numbSim,1)); %cell for number of points
    labelModel='Matern I';
    rSub=rMatern;
elseif choiceModel==2
    %area of hard-core disk
    diskMaternArea=pi*rMatern^2;
    %Matern II
    numberCell=mat2cell(pointCountsII,ones(numbSim,1)); %cell for number of points
    %analytic density
    lambdaSub=1/(diskMaternArea)*(1-exp(-lambda*diskMaternArea)); %Matern II
    labelModel='Matern II';
    rSub=rMatern;
elseif choiceModel==3
    %Triangle
    numberCell=mat2cell(pointCountsTri,ones(numbSim,1)); %cell for number of points
    lambdaSub=mean(pointCountsTri/areaSample);%(empirical estimate)
    labelModel='Triangle';
    rSub=rTri;
end
%create point process object (ie structure array)
ppStructSub=struct('x',xxCellSub,'y',yyCellSub,...
    'n',numberCell,'window',windowSample);

%%%END -- Create point process object -- END %%%

disp(strcat('Finished simulating ',{' '},labelModel,' point process.'));

if ~exist('booleNew','var')
    save('Subset.mat','lambda','xx0','yy0','areaSample','rSub',...
        'lambdaSub','ppStructPoisson','indexCellSub','windowSample',...
        'choiceModel','labelModel','booleDisk');
end

%Checking that all realizations have enough points
if any([ppStructSub(:).n]<=1)
    if (choiceModel==1)||(choiceModel==2)
        warning('There are Matern realizations with one or zero points. Try reducing the parameter rMatern.');
    elseif (choiceModel==3)
        warning('There are triangular realizations with one or zero points. Try reducing the parameter rTri.');
    end
end

%%%Empirical densities of generated point processes
if booleCheck==1
    rSub
    %intensity/density values
    lambdaExact=lambda
    %empirical estimate of density
    lambdaEmp=mean(pointCounts/areaSample)
    
    kk=1;
    xxPoisson=ppStructPoisson(kk).x;
    yyPoisson=ppStructPoisson(kk).y;
    xxSub=xxPoisson(indexCellSub{kk});
    yySub=yyPoisson(indexCellSub{kk});
    
    %%%START -- Plot point processses
    markerSizeNumb=50; %marker size of markers colors
    markerSizeRatio=1;
    figure; hold on;
    set(gcf,'DefaultLineLineWidth',1.2);
    %draw underlying Poisson point process
    plot(xxPoisson,yyPoisson,'ko','MarkerSize',markerSizeNumb/4.0);hold on;
    plot(xxSub,yySub,'rx','MarkerSize',markerSizeNumb/6.0);
    legend('Poisson', labelModel);
    axis square;
    if choiceModel==1
        %analytic density
        lambdaSub
        %empirical estimate of density
        lambdaSubEmp=mean(pointCountsI/areaSample)
    elseif choiceModel==2
        %analytic density
        lambdaSub
        %empirical estimate of density
        lambdaSubEmp=mean(pointCountsII/areaSample)
    else
        %empirical estimate of density (analytic one is not known)
        lambdaSubEmp=lambdaSub
    end   
    
end
