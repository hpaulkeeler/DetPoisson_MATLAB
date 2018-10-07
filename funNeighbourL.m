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

function L=funNeighbourL(xx,yy,lambda,choiceKernel,sigma,theta,N,M)
%Check that M is the right value if it exists
if exist('M','var')
    if (M~=N-1)&&(M>0)
        error('M must be equal to N-1 or zero.');
    end
end
theta=theta(:); %theta needs to be columm vector
xx=xx(:); yy=yy(:); %xx/yy need to be column vectors
meanD=1/sqrt(lambda); %rescale constant
xx=xx/meanD;yy=yy/meanD; %rescale distances

alpha=1; %an additional parameter for the Cauchy kernel

%START - Creation of L matrix
sizeL=length(xx); %width/height of L (ie cardinality of state space)

%START -- Create q (ie quality feature/covariage) vector
%zeroth term
thetaFeature=theta(1)*ones(sizeL,1);
if N>0
    %1st to N th terms
    booleMatrix=logical(ones(sizeL)-eye(sizeL)); %removes distances to self
    indexMatrix=repmat((1:sizeL),sizeL,1);
    %find every pair combination
    xxMatrix=repmat(xx,1,sizeL-1);yyMatrix=repmat(yy,1,sizeL-1);
    indexTemp=flipud(reshape(indexMatrix(booleMatrix),sizeL,sizeL-1));
    %calculate nearest distance for all pairs
    [distNear,indexNear]=sort(sqrt((xxMatrix-xxMatrix(indexTemp)).^2 ...
    +(yyMatrix-yyMatrix(indexTemp)).^2),2);    
    thetaMatrix=repmat(theta(2:end),1,sizeL);
    thetaFeature=thetaFeature ...
        +sum((thetaMatrix(1:N,:)').*distNear(:,1:N),2);
    
    %Run for distances between nearest neighbours
    %N+1 to M th terms
    if exist('M','var')
        %find indices of nearest neighbours
        indexBetween=indexTemp(sub2ind([sizeL,sizeL-1], ...
            repmat((1:sizeL)',1,sizeL-1),indexNear));
        %distance between nearest neighbours neighbours
        distBetween=sqrt((xxMatrix(indexBetween(:,1:M)) ...
            -xxMatrix(indexBetween(:,2:M+1))).^2 ...
            +(yyMatrix(indexBetween(:,1:M)) ...
            -yyMatrix(indexBetween(:,2:M+1))).^2);
        %add contribution
        thetaFeature=thetaFeature ...
            +sum((thetaMatrix(N+1:N+M,:)').*distBetween,2);
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
        SMatrix=(1./1)*exp(-(rrDiffSquared)/sigma^2);
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



% %%% CODE GRAVEYARD -- this code works, but it is slower %%%%
%
%
% %START -- Create q vector
% %%input has to be column vectors for knnsearch
% indexNeighborAll=knnsearch([xx,yy],[xx,yy],'K',(N+1));
% thetaFeature=ones(sizeL,1);
% for ii=1:sizeL
%     indexNearTemp=indexNeighborAll(ii,:); %nearest neighbour index
%     xxTemp=xx(indexNearTemp);yyTemp=yy(indexNearTemp);
%     distNearTemp=sqrt((xx(ii)-xxTemp(2:end)).^2+(yy(ii)-yyTemp(2:end)).^2);
%
%     thetaFeature(ii)=(theta(1:N)'*distNearTemp);
%
%     if exist('M','var')&& (M==N-1)&&(M>0)
%             distBetweenTemp=sqrt((xxTemp(2:M+1)-xxTemp(3:M+2)).^2 ...
%                 +(yyTemp(2:M+1)-yyTemp(3:M+2)).^2);
%             thetaFeature(ii)=thetaFeature(ii) ...
%                 +(theta(N+1:end)'*distBetweenTemp);
%     end
% end
% qVector=exp(thetaFeature);
%
%END -- Create q vector

% %START - Create S matrix with Gaussian kernel
% SMatrix=zeros(sizeL,sizeL);
% for ii=1:sizeL
%     xxDiff1=xx-xx(ii);
%     yyDiff1=yy-yy(ii);
%     xxDiff=xxDiff1; yyDiff=yyDiff1;
%     dd=sqrt((xxDiff).^2+(yyDiff).^2);
%     SMatrix(:,ii)=exp(-(dd.^2/sigma^2)); %Gaussian kernel
% end
% L=(qMatrix').*SMatrix.*qMatrix;

%END - Create S matrix

%Old more complicated Bessel kernel

% LBessel=rho*besselj(constBessel,zDiffScaled)./(zDiffScaled).^((nu+2)/2);
% %need to rescale to ensure that diagonal entries are ones.
% LBesel=(2^constBessel)*gamma(constBessel+1)*LBessel;
% LBesel(1:1+size(LBesel,1):end)=rho; %remove the nan from zero division

%
% elseif
%     rrDiff=sqrt(rrDiffSquared);
%     %%Bessel kernel
%     rho=sigma;
%     %Bessel (simplified) kernel
%     SMatrix=sqrt(rho)*besselj(1,2*sqrt(pi*rho)*rrDiff)./sqrt(pi*rrDiff);
%     %need to rescale to ensure that diagonal entries are ones.
%     SMatrix(1:1+size(SMatrix,1):end)=sqrt(rho); %remove the nan from zero division
