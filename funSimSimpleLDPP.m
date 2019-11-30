% function indexDPP=funSimSimpleLDPP(eigenVectL,eigenValL)
% This function simulates a determinantal point process (DPP) provided a
% L matrix. It was used to produce the results in the paper by Blaszczyszyn
% and Keeler[1].
%
% The algorithm exists in Kulesza and Taskar[2]; see Algorithm 1.
% Also see Algorithm 1 in Lavancier, Moller and Rubak[3].
%
% If you use this code in a publication, please cite the paper by
% Blaszczyszyn and Keeler[1].
%
% Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
% Melbourne, 2018
%
% References:
% [1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
% with network learning applications, 2018.
% [2] Kulesza and Taskar, "Determinantal point processes for machine
% learning",Now Publisers, 2012
% [3] Lavancier, Moller and Rubak, "Determinantal point process models and
% statistical inference", Journal of the Royal Statistical Society --
% Series B, 2015.
function indexDPP=funSimSimpleLDPP(eigenVectL,eigenValL)

%START - Sampling/simulating DPP - START
%Retrieve eigenvalues and eigenvectors
eigenVectK = (eigenVectL); 
eigenValL=abs(diag(eigenValL)); %eigenvalues of L
eigenValK = eigenValL./(1+eigenValL); %eigenvalues of K

%Bernoulli trials (ie coin flips) to determine number of points
booleEigen = (rand(length(eigenValK),1) <= eigenValK);

%number of points in the DPP realization
numbPointsDPP=sum(booleEigen);
%retrieve eigenvectors corresponding to successful Bernoulli trials
spaceV = eigenVectK(:,booleEigen); %subspace V
indexDPP = zeros(numbPointsDPP,1); %index for final DPP configuration

%Loop through for all points
for ii=1:numbPointsDPP
    %Compute probabilities for each point i
    Prob_i=sum(spaceV.^2,2); %sum across rows
    Prob_i=Prob_i/sum(Prob_i); %normalize
    
    %Choose a point (from 1 to numbPoints) using (prob mass function) Prob_i
    uRand=rand(1);
    indexCurrent=find(uRand<=cumsum(Prob_i),1);
    indexDPP(ii)=indexCurrent;
    
    if ii<numbPointsDPP
        %Choose a vector to remove
        jj=find(spaceV(indexCurrent,:),1);
        columnVj=spaceV(:,jj); %j-th column of V
        spaceV=spaceV(:,[1:jj-1 jj+1:end]);
        
        %Update matrix V by removing Vj component from the space
        spaceV=spaceV-columnVj.*spaceV(indexCurrent,:)/columnVj(indexCurrent);
        %Orthonormalize using Householder method
        [spaceV,~]=qr(spaceV,0);
    end
end
%END - Sampling/simulating DPP - END
end



