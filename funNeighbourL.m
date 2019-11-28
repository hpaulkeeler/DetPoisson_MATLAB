% function L=funNeighbourL(xx,yy,lambda,choiceKernel,sigma,theta,N,M)
% This function file creates an L(-ensemble-)matrix, as detailed in the
% paper by Blaszczyszyn and Keeler[1](Section IV).
%
% The quality features or covariates (q_x in the above paper) are based on
% the nearest neighbour distances. The similiarirty matrix (S in the paper)
% which creates replusion among the points, can be formed from either
% Gaussian or Cauchy kernel function.
%
% INPUTS:
% xx and yy are the x and y values of the underlying discrete state space,
% which is usually a realization of a point process on a bounded continuous
% 2-D space such as a square or a disk.
%
% lambda is the point intensity/density (ie average number of points per
% unit area) of the point process, which is used to rescale the distances.
%
% choiceKernel is a variable that takes value 1 (for Gaussian) or 2 (for
% Cauchy) to select the kernel function.
%
% sigma is a parameter of kernel function.
%
% theta is a fitting parameters for the quality features/covariates.
%
% N is the number of neighbouring points.
%
% M is the number of distances between neighbour points. M is optional, but
% when used, M must be equal to zero or N-1.
%
% OUTPUTS:
% An L-ensemble kernel matrix for a determinantal point process on a
% discrete space; see [1] for details.
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

function L=funNeighbourL(xx,yy,lambda,choiceKernel,sigma,theta,N,M)
%Check that M is the right value if it exists
if exist('M','var')
    if (M~=N-1)&&(M>0)
        error('M must be equal to N-1 or zero.');
    end
end
%Convert any matrices into vectors and rescale
theta=theta(:); %theta needs to be columm vector
xx=xx(:); yy=yy(:); %xx/yy need to be column vectors
meanD=1/sqrt(lambda); %rescaling constant
xx=xx/meanD;yy=yy/meanD; %rescale distances

alpha=1; %an additional parameter for the Cauchy kernel

%START - Creation of L matrix
sizeL=length(xx); %width/height of L (ie cardinality of state space)

%START -- Create q (ie quality feature/covariage) vector
%zeroth term
thetaFeature=theta(1)*ones(sizeL,1);
if N>0
    %1st to N th terms
    booleMatrix=logical(ones(sizeL)-eye(sizeL));%removes distances to self
    rowMatrix=repmat((1:sizeL),sizeL,1);
    %find every pair combination
    xxMatrix=repmat(xx,1,sizeL-1);
    yyMatrix=repmat(yy,1,sizeL-1);
    
    %find rows for accessing x/y values
    rowTemp=flipud(reshape(rowMatrix(booleMatrix),sizeL,sizeL-1));
    %access x/y values and reshape vectors
    xxTemp=xxMatrix(rowTemp);
    yyTemp=yyMatrix(rowTemp);
    
    %differences in distances for all pairs
    xxNearDiff=xxMatrix-xxTemp;yyNearDiff=yyMatrix-yyTemp;
    %calculate nearest distance for all pairs
    distNearTemp=hypot(xxNearDiff,yyNearDiff);
    
    %sort distances
    [distNear,indexNear]=sort(distNearTemp,2);
    
    %replicate parameter vector
    thetaMatrix=repmat(theta(2:end),1,sizeL);
    %dot product of parameters and features
    theta_distNear=(thetaMatrix(1:N,:)').*distNear(:,1:N);
    %add contribution fromm parameters and features
    thetaFeature=thetaFeature+sum(theta_distNear,2);
    
    %Run for distances between nearest neighbours
    %N+1 to M th terms
    if exist('M','var')
        %find indices of nearest neighbours
        indexBetweenTemp=sub2ind([sizeL,sizeL-1], ...
            repmat((1:sizeL)',1,sizeL-1),indexNear);
        rowBetween=rowTemp(indexBetweenTemp);
        
        %access x/y values and reshape vectors
        xxTemp=xxMatrix(rowBetween);
        yyTemp=yyMatrix(rowBetween);
        
        %differences in distances for all pairs
        xxNearDiff=xxTemp(:,1:M)-xxTemp(:,2:M+1);
        yyNearDiff=yyTemp(:,1:M)-yyTemp(:,2:M+1);
        %distance between nearest neighbours neighbours
        distBetween=hypot(xxNearDiff,yyNearDiff);
        
        %dot product of parameters and features
        theta_distBetween=(thetaMatrix(N+1:N+M,:)').*distBetween;
        %add contribution from parameters and features
        thetaFeature=thetaFeature+sum(theta_distBetween,2);
    end
    
end
qVector=exp(thetaFeature); %find q vector (ie feature/covariate values)
%END -- Create q vector

%START - Create similarity matrix S
if sigma~=0
    %all squared distances of x/y difference pairs
    xxDiff=bsxfun(@minus,xx,xx'); yyDiff=bsxfun(@minus,yy,yy');
    rrDiffSquared=(xxDiff.^2+yyDiff.^2);
    if choiceKernel==1
        %%Gaussian kernel
        SMatrix=exp(-(rrDiffSquared)/sigma^2);
    else
        
        %%Cauchy kernel
        SMatrix=1./(1+rrDiffSquared/sigma^2).^(alpha+1/2);
    end
else
    SMatrix=eye(sizeL);
end
%END - Create similarity matrix S with Gaussian kernel

%START Create L matrix
qMatrix=repmat(qVector',size(SMatrix,1),1); %q diagonal matrix
L=(qMatrix').*SMatrix.*qMatrix;
%END Create L matrix

end